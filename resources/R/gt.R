box::use(reshape2[melt])
box::use(stringr[str_extract])
box::use(./plot)

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

correct_covmat_for_n <- function(pmat, n){
  pdiff <- pmat[, -1] - pmat[, -ncol(pmat)]

  covmat = cov(pdiff, use = "pairwise.complete.obs")

  for(i in 1:(ncol(covmat)-1)){
    corrected_cov <- mean(pdiff[,i] * pdiff[,i+1], na.rm = T) + mean(pmat[,i+1]*(1-pmat[,i+1])/(n[i+1] - 1), na.rm = T)

    covmat[i,i+1] <- corrected_cov
    covmat[i+1,i] <- corrected_cov
  }

  for(i in 1:ncol(covmat)){
    corrected_var <- mean((pdiff[,i])^2, na.rm = T) - mean(pmat[,i]*(1-pmat[,i])/(n[i] - 1), na.rm = T) - mean(pmat[,i+1]*(1-pmat[,i+1])/(n[i+1] - 1), na.rm = T)
    if(corrected_var < 0) {
      print("Sample size correction makes variance negative. Setting variance to zero")
      covmat[i,i] <- 0
    } else{
      covmat[i,i] <- corrected_var
    }
  }
  return(covmat)
}

gt_n_adjust <- function(pmat, n) {

  if(ncol(pmat) != length(n)){
    stop("Error: Must give a sample size for every timepoint.")
  }

  pdiff <- pmat[, -1] - pmat[, -ncol(pmat)]

  covmat <- cov(pdiff, use = "pairwise.complete.obs")

  covmat <- correct_covmat_for_n(covmat, pmat, n)

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

ssh <- function(x){
  sum(2*x*(1-x))
}

estim_linked_selection_params <- function(pmat, n, mean_ld, weight){

  covmat <- correct_covmat_for_n(pmat, n)

  X <- melt(covmat)
  X <- X[!duplicated(X$value),]

  X$t <- gsub("_", "", str_extract(X$Var1, "_.*_"))
  X$s <- gsub("_", "", str_extract(X$Var2, "_.*_"))

  X$b <- as.numeric(X$t == X$s)

  ssh_vec <- apply(pmat, MARGIN = 2, FUN = ssh)

  ssh_vec <- data.frame(
    t = names(ssh_vec),
    ssh = ssh_vec
  )

  ssh_vec$t <- gsub("_", "", str_extract(ssh_vec$t, "_.*_"))

  X <- merge(X, ssh_vec, by = "t")

  X <- merge(X, ssh_vec, by.x = "s", by.y = "t", suffixes = c("_s", "_t"))

  ssh_1 <- ssh(pmat[,1])

  X$ssh_ratio <- X$ssh_t/ssh_1

  X$a <- X$ssh_ratio*0.5*mean_ld*weight

  return(X)
}
