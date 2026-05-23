box::use(reshape2[melt],
         stringr[str_extract],
         ./plot,
         stats[cov],
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
covmat_from_pmat <- function(pmat, n, correct_for_n = T){
  
  # allele frequency changes between adjacent generations
  pdiff <- freq_increments(pmat)

  # covariance matrices
  covmat <- stats::cov(pdiff, use = "pairwise.complete.obs")

  # correct for sample size, if needed
  if(correct_for_n){
      if(any(n < 2)){
        stop("Error: n must be >= 2")
      }
      
      if(length(n) != (ncol(covmat) + 1)){
        stop("Should be a sample size for every time point.")
      }
    
      # correct overlapping covariances
      for(i in 1:(ncol(covmat)-1)){
        corrected_cov <- mean(pdiff[,i] * pdiff[,i+1], na.rm = T) + mean(pmat[,i+1]*(1-pmat[,i+1])/(n[i+1] - 1), na.rm = T)
        covmat[i,i+1] <- corrected_cov
        covmat[i+1,i] <- corrected_cov
      }
      
      # correct variances
      for(i in 1:ncol(covmat)){
        corrected_var <- mean((pdiff[,i])^2, na.rm = T) - mean(pmat[,i]*(1-pmat[,i])/(n[i] - 1), na.rm = T) - mean(pmat[,i+1]*(1-pmat[,i+1])/(n[i+1] - 1), na.rm = T)
        if(corrected_var < 0) {
          print("Sample size correction makes variance negative. Setting variance to zero")
          covmat[i,i] <- 0
        } else{
          covmat[i,i] <- corrected_var
        }
      }
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

fc <- function(p0, pt) {
  stopifnot(length(p0) == length(pt),
            all(p0 >= 0),
            all(pt >= 0),
            all(p0 <= 1),
            all(pt <= 1),
            all(!is.na(p0)),
            all(!is.na(pt)))
  
  L <- length(p0)
  
  q0 <- 1 - p0
  
  qt <- 1 - pt
  
  fsum <- ((p0 - pt) ^ 2) / ((p0 + pt) / 2 - p0 * pt) + ((q0 - qt) ^ 2) /
    ((q0 + qt) / 2 - q0 * qt)
  
  #return(fsum)
  return(sum(fsum, na.rm = T) / sum(!is.nan(fsum)))
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
waples_ne <- function(s_low, s_high, Fc, t, S0, St){
  r <- c(ncne(s_low), ncne(s_high))
  ne <- (r*t - 2)/(2*r*(Fc - 1/(2*S0) - 1/(2*St)))
  return(ne)
}
