box::use(reshape2[melt],
         stringr[str_extract],
         ./plot,
         stats[...],
         dplyr[...]
         )

#' Calculate linked selection ratios from covariance matrix
#'
#' @param covmat 
#'
#' @return
#' @export
#'
#' @examples
gt_from_covmat <- function(covmat) {
  covmat[lower.tri(covmat)] <- NA

  covariances <- covmat[upper.tri(covmat)]
  variances <- diag(covmat)

  sums_var <- cumsum(variances)

  sums_cov = c(0)
  sums_abs_cov = c(0)
  sums_pos_cov = c(0)
  sums_neg_cov = c(0)
  for(j in 2:nrow(covmat)){
    cov_sub <- covariances[1:(((j-1)*(j))/2)]

    sum_abs_cov <- 2*sum(abs(cov_sub), na.rm = T)
    sum_cov <- 2*sum(cov_sub, na.rm = T)
    sum_neg_cov <- 2*sum(abs(cov_sub[(cov_sub < 0)]), na.rm = T)
    sum_pos_cov <- 2*sum(cov_sub[(cov_sub > 0)], na.rm = T)

    sums_cov = c(sums_cov, sum_cov)
    sums_abs_cov = c(sums_abs_cov, sum_abs_cov)
    sums_pos_cov = c(sums_pos_cov, sum_pos_cov)
    sums_neg_cov = c(sums_neg_cov, sum_neg_cov)
  }

  result = data.frame(
    gen = 1:nrow(covmat),
    gt = sums_cov/(sums_var + sums_cov),
    gt_abs = sums_abs_cov/(sums_var + sums_abs_cov),
    gt_abs_2 = sums_abs_cov/(sums_var + sums_cov),
    gt_pos = sums_pos_cov/(sums_var + sums_pos_cov),
    sum_cov = sums_cov,
    sum_abs_cov = sums_abs_cov,
    sum_pos_cov = sums_pos_cov,
    sum_neg_cov = sums_neg_cov,
    sum_var = sums_var,
    ratio = sums_pos_cov/(sums_pos_cov + sums_neg_cov)
  )
  return(result)
}

gt_from_pmat <- function(dps, calcabs = FALSE, output_name) {
  results = NULL
  col_names <- colnames(dps[[1]])
  for (i in 1:ncol(dps[[1]])) {
    print(i)
    dp = NULL
    for (k in 1:length(dps)) {
      print(k)
      dp <- cbind(dp, (dps[[k]])[,i])
    }

    print("Calculating covariance...")
    covmat <- cov(dp)
    covmat[lower.tri(covmat)] <- NA

    plot$plot_var_cov_matrix(covmat, paste(output_name, "_", i, ".pdf", sep = ""))

    print("Extracting covariances and variances...")
    covariances <- covmat[upper.tri(covmat)]
    variances <- diag(covmat)

    print("Cumulative sum of variances...")
    sums_var <- cumsum(variances)

    print("Cumulative sum of covariances...")
    sums_cov = c(0)
    sums_abs_cov = c(0)
    sums_pos_cov = c(0)
    sums_neg_cov = c(0)
    for (j in 2:nrow(covmat)) {
      print(j)
      cov_sub <- covariances[1:(((j - 1) * (j)) / 2)]

      sum_abs_cov <- 2 * sum(abs(cov_sub), na.rm = T)
      sum_cov <- 2 * sum(cov_sub, na.rm = T)
      sum_neg_cov <- 2 * sum(abs(cov_sub[(cov_sub < 0)]), na.rm = T)
      sum_pos_cov <- 2 * sum(cov_sub[(cov_sub > 0)], na.rm = T)

      sums_cov = c(sums_cov, sum_cov)
      sums_abs_cov = c(sums_abs_cov, sum_abs_cov)
      sums_pos_cov = c(sums_pos_cov, sum_pos_cov)
      sums_neg_cov = c(sums_neg_cov, sum_neg_cov)
    }

    print("Calculating gt...")
    if (calcabs == TRUE) {
      gts = sums_abs_cov / (sums_var + sums_abs_cov)
    } else {
      gts = sums_cov / (sums_var + sums_cov)
    }

    print("Saving result...")
    result = data.frame(
      gen = 1:nrow(covmat),
      gt = gts,
      sum_cov = sums_cov,
      sum_abs_cov = sums_abs_cov,
      sum_pos_cov = sums_pos_cov,
      sum_neg_cov = sums_neg_cov,
      ratio = sums_neg_cov / (sums_pos_cov + sums_neg_cov),
      sum_var = sums_var,
      replicate = i
    )

    results = rbind(results, result)
  }
  return(results)
}

#' Allele frequency changes between adjacent generations
#'
#' @param pmat matrix of allele frequencies, timepoints ascending from left column to right column
#'
#' @return
#' @export
#'
#' @examples
freq_increments <- function(pmat){
  pmat[, -1] - pmat[, -ncol(pmat)]
}


#' Standardize covmat by 
#'
#' @param pmat excluding the last time point
#' @param covmat 
#'
#' @return
#' @export
#'
#' @examples
standard_cov_by_het <- function(pmat, covmat) {
  stopifnot(ncol(pmat) == ncol(covmat))
  half_het_sums <- 0.5*apply(pmat, MARGIN = 2, FUN = sum_of_het)
  stopifnot(all(half_het_sums > 0))
  # standaridize variances
  diag(covmat) <- diag(covmat)/half_het_sums
  # standardize covariances
  for(i in 1:nrow(covmat)){
    for(j in 1:ncol(covmat)){
      if(i != j){
        covmat[i,j] <- covmat[i,j]/half_het_sums[min(c(i,j))]
      }
    }
  }
  return(covmat)
}

#' Calculate covariances from allele frequency matrix
#'
#' @param pmat matrix of allele frequencies, first column is initial generation
#' @param n optional, vector of sample size, number of chromosomes by site
#' @param correct_for_n boolean, whether to correct for sample size
#'
#' @return
#' @export
#'
#' @examples
covmat_from_pmat <- function(pmat, n = NULL, correct_for_n = T, standard_by_het = F, input_asin_trans = F, procedure = "none", windows = NULL){
  # allele frequency changes between adjacent generations
  pdiff <- freq_increments(pmat)
  # randomly swap signs of allele frequency changes
  pdiff <- sign_permute_increments(pdiff, procedure = procedure, windows = windows)
  # ensure pdiff is a matrix even with a single increment (2-column pmat)
  if(is.null(dim(pdiff))){
    pdiff <- matrix(pdiff, ncol = 1)
  }
  # covariance matrices
  covmat <- stats::cov(pdiff, use = "pairwise.complete.obs")
  # correct for sample size, if needed
  if(correct_for_n){
    if(any(n < 2)){
      stop("All n must be >= 2")
    }
    if(length(n) != (ncol(covmat) + 1)){
      stop("Should be a sample size for every time point.")
    }
    ### ARCSIN SQRT TRANSFORMED FREQUENCIES ###
    if(input_asin_trans){
      # check transformed data are actually input
      if(all(pmat[,1] < pi)){
        warning("Are you sure values are arcsin transformed?")
      }
      # correct overlapping covariances
      # seq_len() yields integer(0) when ncol(covmat) == 1, so this is skipped
      for(i in seq_len(ncol(covmat) - 1)){
        corrected_cov <- covmat[i,i+1] + (1/n[i+1])
        covmat[i,i+1] <- corrected_cov
        covmat[i+1,i] <- corrected_cov
      }
      # correct variances
      for(i in seq_len(ncol(covmat))){
        corrected_var <- covmat[i,i] - (1/n[i]) - (1/n[i+1])
        if(corrected_var < 0) {
          print("Sample size correction makes variance negative. Setting variance to zero")
          covmat[i,i] <- 0
        } else{
          covmat[i,i] <- corrected_var
        }
      }
      ###RAW FREQUENCIES###
    } else {
      # correct overlapping covariances
      # seq_len() yields integer(0) when ncol(covmat) == 1, so this is skipped
      for(i in seq_len(ncol(covmat) - 1)){
        corrected_cov <- mean(pdiff[,i] * pdiff[,i+1], na.rm = T) + mean(pmat[,i+1]*(1-pmat[,i+1])/(n[i+1] - 1), na.rm = T)
        covmat[i,i+1] <- corrected_cov
        covmat[i+1,i] <- corrected_cov
      }
      # correct variances
      for(i in seq_len(ncol(covmat))){
        corrected_var <- mean((pdiff[,i])^2, na.rm = T) - mean(pmat[,i]*(1-pmat[,i])/(n[i] - 1), na.rm = T) - mean(pmat[,i+1]*(1-pmat[,i+1])/(n[i+1] - 1), na.rm = T)
        if( corrected_var < 0 | is.na(corrected_var) ) {
          print("Sample size correction makes variance negative. Setting variance to zero")
          covmat[i,i] <- 0
        } else{
          covmat[i,i] <- corrected_var
        }
      }
    }
  }
  
  if(standard_by_het){
    covmat <- standard_cov_by_het(pmat[,-ncol(pmat), drop = FALSE], covmat)
  }
  
  return(covmat)
}

#' Calculate covariances across a pair of environments
#'
#' @param pmat1 matrix of allele frequencies in population 1, rows = variants, columns = replicates
#' @param pmat2 matrix of allele frequencies in population 2, rows = variants, columns = replicates
#' @param pop1 
#' @param pop2 
#' @param gen0 
#'
#' @return
#' @export
#'
#' @examples
covmat_pop_pair <- function(pmat1, pmat2, pop1, pop2, gen0){
  # allele frequency changes between adjacent generations
  pdiff1 <- pmat1 - gen0
  pdiff2 <- pmat2 - gen0
  
  covmat <- cov(cbind(pdiff1, pdiff2), use = "pairwise.complete.obs")
  
  # mean product of std of between population pairs
  rep_var <- diag(covmat)
  rep_var_1 <- rep_var[1:ncol(pmat1)]
  rep_var_2 <- rep_var[(ncol(pmat1) + 1):length(rep_var)]
  mean_pair_std <- mean(sqrt(rep_var_1 %o% rep_var_2))
  
  bw_pop_cov <- covmat[grep(pop1, rownames(covmat)), grep(pop2, colnames(covmat))]
  
  return(mean(bw_pop_cov)/mean_pair_std)
}


#' Pairwise geometric mean
#'
#' @param my_vector 
#'
#' @return
#' @export
#'
#' @examples
geom_pairwise_mean <- function(my_vector) {
  stopifnot(length(my_vector)>1, is.numeric(my_vector))
  n <- length(my_vector)
  var_sum <- 0
  for (i in 1:(n - 1)) {
    var_a <- my_vector[i]
    for (j in (i + 1):n) {
      var_b <- my_vector[j]
      var_sum <- var_sum + sqrt(var_a*var_b)
    }
  }
  var_mean <- var_sum/((n^2 - n)/2)
  return(var_mean)
}

#' Calculate convergence correlation within environment
#'
#' @param pdiff frequency increments for all replicates of one site, one time point
#'
#' @return
#' @export
#'
#' @examples
conv_cor_wn_env <- function(pdiff){
  stopifnot(ncol(pdiff) > 1, nrow(pdiff) > 1)
  covmat <- cov(pdiff, use = "pairwise.complete.obs")
  rep_var <- diag(covmat)
  rep_cov <- covmat[upper.tri(covmat)]
  numerator <- sum(2*rep_cov, na.rm = T)/(2*length(rep_cov))
  denominator <- geom_pairwise_mean(rep_var)
  stopifnot(denominator > 0, abs(numerator) <= denominator)
  return(c(numerator,denominator))
}


#' Replicate G(t): proportion of variance in dp due to shared selection pressure
#'
#' @param pdiff 
#' @param rep_labels 
#' @param time_labels 
#'
#' @return
#' @export
#'
#' @examples
replicate_gt <- function(pdiff, rep_labels, time_labels){
  
  stopifnot(ncol(pdiff) > 1, 
            nrow(pdiff) > 1, 
            length(rep_labels) == ncol(pdiff), 
            length(time_labels) == ncol(pdiff),
            length(unique(rep_labels)) >= 2,
            length(unique(time_labels)) >= 2)
  
  rep_eq_bool_mat <- outer(rep_labels, rep_labels, "==")
  
  time_eq_bool_mat <- outer(time_labels, time_labels, "==")
  
  covmat <- cov(pdiff, use = "pairwise.complete.obs")
  
  # total variance within replicates
  var_only_mat <- covmat[rep_eq_bool_mat]
  total_var <- sum(var_only_mat)
  total_var_pos <- sum(var_only_mat[(var_only_mat > 0)])
  total_var_abs <- sum(abs(var_only_mat))
  var_trace <- covmat[(rep_eq_bool_mat & time_eq_bool_mat)]
  total_var_trace <- sum(var_trace)
  n_var <- length(unique(rep_labels))
  
  # total covariance across different replicates AND time points
  cov_only_mat <- covmat[(!rep_eq_bool_mat) & (!time_eq_bool_mat)]
  total_covar <- sum(cov_only_mat)
  total_covar_abs <- sum(abs(cov_only_mat))
  total_covar_pos <- sum(cov_only_mat[(cov_only_mat > 0)])
  n_covar <- (n_var^2 - n_var)/2
  
  mean_var = total_var/n_var
  mean_var_pos = total_var_pos/n_var
  mean_var_abs = total_var_abs/n_var
  mean_var_trace = total_var_trace/n_var
  mean_covar = total_covar/n_covar
  mean_covar_pos = total_covar_pos/n_covar
  mean_covar_abs = total_covar_abs/n_covar
  
  stopifnot(total_covar <= total_covar_pos,
            mean_covar <= mean_covar_pos,
            total_var_pos <= total_var_abs,
            total_covar < sum(abs(covmat)), 
            total_var < sum(abs(covmat)),
            total_covar <= total_covar_abs,
            all.equal(sum(var_only_mat)+sum(covmat[!rep_eq_bool_mat]), 
                      sum(covmat)),
            mean_var >= 0, 
            total_var >= 0,
            total_var_trace >= 0,
            mean_var_trace >= 0)
  
  return(data.frame(
    mean_var = mean_var,
    mean_var_pos = mean_var_pos,
    mean_var_abs = mean_var_abs,
    mean_var_trace = mean_var_trace,
    mean_covar = mean_covar,
    mean_covar_pos = mean_covar_pos,
    mean_covar_abs = mean_covar_abs,
    total_var = total_var,
    total_var_pos = total_var_pos,
    total_var_abs = total_var_abs,
    total_var_trace = total_var_trace,
    total_covar = total_covar,
    total_covar_pos = total_covar_pos,
    total_covar_abs = total_covar_abs,
    n_var = n_var,
    n_covar = n_covar
  ))

}


#' Rolling sum of matrix, starting from top left, sum matrix elements in nested squares
#'
#' @param mat 
#'
#' @return 
#' @export
#'
#' @examples
rolling_matrix_sum <- function(mat){
  stopifnot(nrow(mat) == ncol(mat))
  rolling_sum <- c()
  for(i in 1:ncol(mat)){
    rolling_sum <- c(rolling_sum, sum(mat[1:i,1:i]))
  }
  stopifnot(length(rolling_sum) == ncol(mat))
  return(rolling_sum)
}

#' G statistic, corrected for selection increasing variance in allele frequency change
#'
#' @param t 
#' @param p0 
#' @param N 
#' @param pmat 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
g_prime <- function(times, pmat, N, n, take_abs = F){
  stopifnot(all(times >= 0), N > 2)
  covmat <- covmat_from_pmat(pmat, n)
  ep0 <- mean(pmat[,1]*(1-pmat[,1]), na.rm = T)
  if(take_abs){
    var_roll <- rolling_matrix_sum(abs(covmat))
  } else {
    var_roll <- rolling_matrix_sum(covmat)
  }
  stopifnot(length(times) == length(var_roll))
  return(1 - (times*ep0)/(2*N*var_roll))
}



#' Sign permute a block of allele frequencies
#'
#' @param pdiff 
#' @param procedure 
#' @param windows 
#'
#' @returns
#' @export
#'
#' @examples
sign_permute_increments <- function(pdiff, procedure = "none", windows = NULL){
  stopifnot(ncol(pdiff)>=2, nrow(pdiff)>=2)
  if(procedure == "window"){
    # within a window, randomly multiply allele frequency change by +1 or -1
    pdiff$window <- windows
    pdiff <- pdiff %>%
      group_by(window) %>%
      mutate(across(where(is.numeric), ~ .x * sample(c(1, -1), size = 1))) %>%
      ungroup() %>%
      dplyr::select(-window) %>%
      as.data.frame()
  } else if (procedure == "cell"){
    # randomly multiply every cell frequency change by +1 or -1
    pdiff <- pdiff %>%
      mutate(across(where(is.numeric), ~ .x * sample(c(1, -1), n(), replace = TRUE)))
  } else if (procedure == "genome"){
    # randomly multiply entire columns by +1 or -1
    sign_permutations <- sample(c(-1,1), size = ncol(pdiff), replace = T)
    # Create the diagonal matrix from vector
    P <- diag(sign_permutations)
    # Perform matrix multiplication to flip all increments at once
    pdiff <- as.matrix(pdiff) %*% P
  } else if (procedure == "none") {
    # no flipping
  }
  else {
    stop("Procedure not valid. Choose from list of options.")
  }
  return(pdiff)
}

gt_n_adjust <- function(pmat, n) {

  if(ncol(pmat) != length(n)){
    stop("Error: Must give a sample size for every timepoint.")
  }

  pdiff <- pmat[, -1] - pmat[, -ncol(pmat)]

  covmat <- cov(pdiff, use = "pairwise.complete.obs")

  covmat <- covmat_from_pmat(covmat, pmat, n)

  covmat[lower.tri(covmat)] <- NA
  covariances <- covmat[upper.tri(covmat)]
  variances <- diag(covmat)

  sums_var <- cumsum(variances)

  sums_cov = c(0)
  sums_abs_cov = c(0)
  sums_pos_cov = c(0)
  sums_neg_cov = c(0)
  for(j in 2:nrow(covmat)){
    cov_sub <- covariances[1:(((j-1)*(j))/2)]

    sum_abs_cov <- 2*sum(abs(cov_sub), na.rm = T)
    sum_cov <- 2*sum(cov_sub, na.rm = T)
    sum_neg_cov <- 2*sum(abs(cov_sub[(cov_sub < 0)]), na.rm = T)
    sum_pos_cov <- 2*sum(cov_sub[(cov_sub > 0)], na.rm = T)

    sums_cov = c(sums_cov, sum_cov)
    sums_abs_cov = c(sums_abs_cov, sum_abs_cov)
    sums_pos_cov = c(sums_pos_cov, sum_pos_cov)
    sums_neg_cov = c(sums_neg_cov, sum_neg_cov)
  }

  result = data.frame(
    gen = 1:nrow(covmat),
    gt = sums_cov/(sums_var + sums_cov),
    gt_abs = sums_abs_cov/(sums_var + sums_abs_cov),
    sum_cov = sums_cov,
    sum_abs_cov = sums_abs_cov,
    sum_pos_cov = sums_pos_cov,
    sum_neg_cov = sums_neg_cov,
    sum_var = sums_var,
    ratio = sums_pos_cov/(sums_pos_cov + sums_neg_cov)
  )

  return(result)
}

#' Sum of heterozygosity 
#'
#' @param x vector of allele frequencies, one per locus
#'
#' @return single number
#' @export
#'
#' @examples
sum_of_het <- function(x){
  sum(2*x*(1-x), na.rm = T)
}


#' Linear interpolation of recombination map 
#'
#' @param x0 vector of positions in bp of target loci
#' @param x1 position in bp, left ref locus
#' @param y1 postion in cM, left ref locus
#' @param x2 position in bp, right ref locus
#' @param y2 position in cM, right ref locus
#'
#' @return vector of position in cM of target loci
#' @export
#'
#' @examples
linear_interpolation <- function(x0, x1, y1, x2, y2) {
  
  # check input format
  if(any(c(x1,y1,x2,y2) < 0)){
    stop("Ref points should all be positive numbers.")
  }
  
  if(x2 < x1){
    stop("x2 should be > x1")
  }
  
  if(y2 < y1){
    stop("y2 should be > y1")
  }
  
  if(any(x0 < x1) | any(x0 > x2)){
    stop("x0 is outsite of [x1, x2]")
  }
  
  # slope
  m = (y2 - y1)/(x2 - x1)
  
  # solve for intercept
  b = y1 - m*x1
  
  # now interpolate
  y0 = m*x0 + b
  
  return(y0)
}


#' sum of heterozygosity across sites, 
#'
#' @param pmat 
#'
#' @return
#' @export
#'
#' @examples
sum_of_het_by_t <- function(pmat){
  
  sum_of_het_vec <- apply(pmat, MARGIN = 2, FUN = sum_of_het)
  
  sum_of_het_vec <- data.frame(
    t = names(sum_of_het_vec),
    sum_of_het = sum_of_het_vec
  )
  
  sum_of_het_vec$t <- gsub("_", "", str_extract(sum_of_het_vec$t, "_.*_"))
  
  return(sum_of_het_vec)
}

#' Compile table of values in windows for Buffalo and Coop 2019 Ne and VA(1) estimation
#'
#' @param pmat 
#' @param n 
#' @param mean_ld 
#' @param weight 
#' @param sum_of_het_vec 
#' @param sum_of_het_1 
#'
#' @return
#' @export
#'
#' @examples
estim_linked_selection_params <- function(pmat, sum_of_het_vec, sum_of_het_1, n, mean_ld, weight){

  covmat <- covmat_from_pmat(pmat, n)

  X <- melt(covmat)
  X <- X[!duplicated(X$value),]

  X$t <- gsub("_", "", str_extract(X$Var2, "_.*_"))
  X$s <- gsub("_", "", str_extract(X$Var1, "_.*_"))

  X$b <- as.numeric(X$t == X$s)

  X <- merge(X, sum_of_het_vec, by = "t")

  X <- merge(X, sum_of_het_vec, by.x = "s", by.y = "t", suffixes = c("_t", "_s"))

  #sum_of_het_1 <- sum_of_het(pmat[,1])

  X$sum_of_het_ratio <- X$sum_of_het_s/sum_of_het_1
  
  #X$sum_of_het_ratio <- X$sum_of_het_s/X$sum_of_het_t

  X$a <- X$sum_of_het_ratio*0.5*mean_ld*weight
  
  X$mean_ld <- mean_ld
  
  X$weight <- weight
  
  X$sum_of_het_1 <- sum_of_het_1

  return(X)
}



#' Title
#'
#' @param pmat 
#' @param sum_of_het_vec 
#' @param sum_of_het_1 
#' @param n 
#' @param mean_ld 
#' @param weight 
#'
#' @return
#' @export
#'
#' @examples
estim_linked_selection_params_new <- function(pmat, n, sum_of_het_vec, sum_of_het_1, mean_ld, weight, correct_for_n = T){
  
  covmat <- covmat_from_pmat(pmat, n = n, correct_for_n = correct_for_n, standard_by_het = T)
  
  if(is.null(rownames(covmat))){
    rownames(covmat) = colnames(pmat)[2]
    colnames(covmat) = colnames(pmat)[2]
  }
  
  X <- reshape2::melt(covmat)
  X <- X[!duplicated(X$value),]
  
  X$t <- gsub("_", "", str_extract(X$Var2, "_.*_"))
  X$s <- gsub("_", "", str_extract(X$Var1, "_.*_"))
  
  X$b <- as.numeric(X$t == X$s)
  
  X <- merge(X, sum_of_het_vec, by = "t")
  
  X <- merge(X, sum_of_het_vec, by.x = "s", by.y = "t", suffixes = c("_t", "_s"))
  
  #sum_of_het_1 <- sum_of_het(pmat[,1])
  
  X$sum_of_het_ratio <- X$sum_of_het_s/sum_of_het_1
  
  #X$sum_of_het_ratio <- X$sum_of_het_1/X$sum_of_het_s
  
  #X$sum_of_het_ratio <- X$sum_of_het_s/X$sum_of_het_t
  
  #X$sum_of_het_ratio <- X$sum_of_het_t/sum_of_het_1
  
  X$a <- (X$sum_of_het_ratio)*0.5*mean_ld*weight
  
  X$mean_ld <- mean_ld
  
  X$weight <- weight
  
  X$sum_of_het_1 <- sum_of_het_1
  
  return(X)
}

#' Allele frequency trajectory simulation, THE SIGNATURE OF POSITIVE SELECTION ON STANDING GENETIC VARIATION (PRZEWORSKI 2005)
#'
#' @param p0 
#' @param N 
#' @param s 
#' @param t 
#'
#' @return
#' @export
#'
#' @examples
simulate_freq_traj <- function(p0, N, s, t){
 traj <- c(p0)
 p_i <- p0
 dt <- 1/(4*N)
 for(i in 1:t){
   rng <- runif(1)
   if(rng > 0.5){
     p_next <- p_i + 2*N*s*p_i*(1 - p_i)*dt - sqrt(p_i*(1-p_i)*dt)
   } else{
     p_next <- p_i + 2*N*s*p_i*(1 - p_i)*dt + sqrt(p_i*(1-p_i)*dt)
   }
   traj <- c(traj, p_next)
   p_i <- p_next
 }
 return(traj)
}

#' Simulate allele frequency estimation, Buffalo and Coop 2020 PNAS
#'
#' @param p0 
#' @param n 
#' @param d 
#'
#' @return
#' @export
#'
#' @examples
simulate_af_est <- function(p, n, d){
  d_obs <- rpois(1, d)
  if (d_obs == 0) return(NA_real_)
  x_obs <- rbinom(1, size = 2*n, prob = p)
  r_obs <- rbinom(1, size = d_obs, prob = x_obs/(2*n))
  p_obs <- r_obs/d_obs
  return(p_obs)
}

#' Nc/Ne ratio from selfing rate, Pollak 1987
#'
#' @param s 
#'
#' @return
#' @export
#'
#' @examples
ncne <- function(s){
  stopifnot(is.numeric(s), s >= 0, s <= 1)
  1/(1-s/2)
}

#' Estimating standardized variance from allele frequencies
#'
#' @param p0 
#' @param pt 
#'
#' @return
#' @export
#'
#' @examples
fc <- function(p0, pt) {
  stopifnot(length(p0) == length(pt),
            all(p0 >= 0, na.rm = T),
            all(pt >= 0, na.rm = T),
            all(p0 <= 1, na.rm = T),
            all(pt <= 1, na.rm = T))
  
  #L <- length(p0)
  
  fsum_num <- (p0 - pt)^2
  
  #fsum_denom <- ((p0 + pt) / 2 - p0 * pt)
  fsum_denom <- p0*(1-p0)
  return(sum(fsum_num, na.rm = T)/sum(fsum_denom, na.rm = T))
  #return(fsum)
  #return(sum(fsum, na.rm = T) / sum(!is.nan(fsum)))
}

#' Estimate Ne, given estimate of N/Ne from selfing rate, Waples 1989
#'
#' @param s_low 
#' @param s_high 
#' @param Fc 
#' @param t 
#' @param S0 
#' @param St 
#'
#' @return
#' @export
#'
#' @examples
waples_ne <- function(r_low, r_high, p0, pt, S0, St, t){
  #r <- c(ncne(s_low), ncne(s_high))
  r <- c(r_low, r_high)
  Fc <- fc(p0, pt)
  ne <- (r*t - 2)/(2*r*(Fc - 1/(S0) - 1/(St)))
  return(ne)
}


#' Function to calculate LD decay distances, Hill and Weir 1988
#'
#' @param d distance between two loci in bp
#' @param n sample size
#' @param C Population-scaled recombination rate: 4Nc  
#'
#' @return
#' @export
#'
#' @examples
hill_weir_r2 <- function(d, n, C) {
  part1 <- (10 + C * d) / ((2 + C * d) * (11 + C * d))
  part2 <- 1 + ((3 + C * d) * (12 + 12 * C * d + (C * d) ^ 2)) / (n * (2 + C *
                                                                         d) * (11 + C * d))
  r2 <- part1 * part2
  return(r2)
}
