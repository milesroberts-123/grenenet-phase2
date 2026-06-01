library(testthat)

# Source the functions under test
box::use(./misc)
box::use(./gt)
box::use(./plot)

describe("replicate_gt", {
  it("denominator of expectation scales with replicates", {
    set.seed(123)
    pdiff <- matrix(runif(30*10, min = -1, max = 1), nrow = 30, ncol = 10)
    rep_labels <- c(1,1,2,2,3,3,4,4,5,5)
    time_labels <- c(1,2,1,2,1,2,1,2,1,2)
    result <- replicate_gt(pdiff, rep_labels, time_labels)
    expect_equal(result$n_covar, 10)
    expect_equal(result$n_var, 5)
  })
})

describe("geom_pairwise_mean", {
  it("rejects vector with 1 element", {
    expect_error(gt$geom_pairwise_mean(c(10)))
  })
  
  it("rejects non-numeric vector", {
    expect_error(gt$geom_pairwise_mean(c("hello", "world", "blake")))
  })
  
  it("returns 1 when input is all 1's", {
    expect_equal(gt$geom_pairwise_mean(rep(1, times = 10)), 1)
  })
})

describe("fitfreq", {
  it("returns 0.5 when q=0.5, h=0.5, s=0", {
    expect_equal(misc$fitfreq(0.5, 0.5, 0), 0.5)
  })

  it("returns q when s=0 regardless of h", {
    expect_equal(misc$fitfreq(0.2, 0.0, 0), 0.2)
    expect_equal(misc$fitfreq(0.8, 1.0, 0), 0.8)
  })

  it("returns q1>q0 when h=0 and s!=0 (dominance coefficient is zero)", {
    expect_gt(misc$fitfreq(0.3, 0, 0.1), 0.3)
  })

  it("throws error or returns valid value for q in [0,1]", {
    expect_equal(misc$fitfreq(0.0, 0.5, 0.1), 0.0)
    expect_equal(misc$fitfreq(1.0, 0.5, 0.1), (1 + 0.1) / (1 + 0.1))
  })
})

describe("freq_increments", {
  it("correctly computes differences between adjacent generations", {
    pmat <- matrix(c(0.1, 0.2, 0.3, 0.4,
                     0.2, 0.3, 0.4, 0.5), nrow = 2, ncol = 4, byrow = T)
    expected <- matrix(c(0.1, 0.1, 0.1,
                         0.1, 0.1, 0.1), nrow = 2, ncol = 3, byrow = T)
    result <- gt$freq_increments(pmat)
    expect_true(all(abs(result  - expected) < 1e-9))
  })

  it("returns a matrix with one fewer column than input", {
    pmat <- matrix(1:6, nrow = 2, ncol = 3)
    result <- gt$freq_increments(pmat)
    expect_equal(ncol(result), ncol(pmat) - 1)
    expect_equal(nrow(result), nrow(pmat))
  })

  it("handles 1-row matrix", {
    pmat <- matrix(c(0.1, 0.3, 0.6), nrow = 1)
    result <- gt$freq_increments(pmat)
    expected <- matrix(c(0.2, 0.3), nrow = 1)
    expect_true(all(abs(result - expected) < 1e-9))
  })
})

describe("sum_of_het", {
  it("correctly computes sum of 2*x*(1-x) for a vector", {
    expect_equal(gt$sum_of_het(c(0.5)), 0.5)
    expect_equal(gt$sum_of_het(c(0.0, 1.0)), 0)
    expect_equal(gt$sum_of_het(c(1/3, 1/3, 1/3)), 2 * (1/3) * (2/3) * 3)
  })

  it("returns zero for fixation states (all 0 or all 1)", {
    expect_equal(gt$sum_of_het(c(0, 0, 0)), 0)
    expect_equal(gt$sum_of_het(c(1, 1, 1)), 0)
  })
})

describe("get_upper_tri", {
  it("sets lower triangle and diagonal to NA", {
    mat <- matrix(1:9, nrow = 3)
    result <- plot$get_upper_tri(mat)
    expect_true(all(is.na(result[lower.tri(result)])))
    expect_true(all(is.na(diag(result))))
  })

  it("preserves upper triangle values", {
    mat <- matrix(c(1, 2, 0, 4, 5, 6, 7, 8, 9), nrow = 3)
    result <- plot$get_upper_tri(mat)
    expect_true(result[1, 2] == 4)
    expect_true(result[1, 3] == 7)
    expect_true(result[2, 3] == 8)
  })
})

describe("gt_from_covmat", {
  it("returns 0 gt for 2x2 matrix with only one off-diagonal", {
    covmat <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)
    result <- gt$gt_from_covmat(covmat)
    expect_equal(nrow(result), 2)
    # First row: sums_cov=0, sums_var=1 => gt=0/(1+0)=0
    expect_equal(result$gt[1], 0)
  })

  it("handles symmetric covariance matrix", {
    covmat <- matrix(c(4, 1, 1, 9), nrow = 2)
    result <- gt$gt_from_covmat(covmat)
    expect_equal(result$sum_var[2], 13)  # 4+9
    expect_equal(result$sum_cov[2], 2)
    expect_equal(result$sum_abs_cov[2], 2)
  })
})

describe("rm_na_after_na", {
  it("returns unchanged vector when no NA present", {
    x <- c(1, 2, 3, 4, 5)
    result <- misc$rm_na_after_na(x)
    expect_true(all(result == x))
  })

  it("na-truncates from first NA position", {
    x <- c(1, 2, NA, 4, 5)
    result <- misc$rm_na_after_na(x)
    expected <- c(1, 2, NA, NA, NA)
    expect_equal(result, expected)
  })

  it("handles vector with NA at start", {
    x <- c(NA, 2, 3)
    result <- misc$rm_na_after_na(x)
    expect_true(all(is.na(result)))
  })

  it("returns vector unchanged when NA only at the end", {
    x <- c(1, 2, NA)
    result <- misc$rm_na_after_na(x)
    expected <- c(1, 2, NA)
    expect_equal(result, expected)
  })
})

describe("sum_of_het_by_t", {
  it("returns data frame with columns t and sum_of_het", {
    pmat <- matrix(c(rep(0.5, 15)), nrow = 5, ncol = 3,
                     dimnames = list(NULL, c("gen_1", "gen_2", "gen_3")))
    result <- gt$sum_of_het_by_t(pmat)
    expect_true(is.data.frame(result))
    expect_true("t" %in% names(result))
    expect_true("sum_of_het" %in% names(result))
    expect_equal(nrow(result), 3)
  })

  it("computes correct sum_of_het values", {
    pmat <- matrix(c(rep(0.5, 15)), nrow = 5, ncol = 3,
                     dimnames = list(NULL, c("gen_1", "gen_2", "gen_3")))
    result <- gt$sum_of_het_by_t(pmat)
    # 5 rows * 2*0.5*(1-0.5) = 5 * 0.5 = 2.5
    expect_equal(result$sum_of_het, rep(2.5, 3))
  })

  it("computes 0 for fixation states (all 0 or all 1)", {
    pmat_fix0 <- matrix(rep(0, 12), nrow = 4, ncol = 3,
                        dimnames = list(NULL, c("gen_1", "gen_2", "gen_3")))
    pmat_fix1 <- matrix(rep(1, 12), nrow = 4, ncol = 3,
                        dimnames = list(NULL, c("gen_1", "gen_2", "gen_3")))
    result0 <- gt$sum_of_het_by_t(pmat_fix0)
    result1 <- gt$sum_of_het_by_t(pmat_fix1)
    expect_equal(result0$sum_of_het, rep(0, 3))
    expect_equal(result1$sum_of_het, rep(0, 3))
  })

  it("produces different values for different allele frequencies", {
    pmat <- matrix(c(rep(0.5, 10),
                     rep(0.1, 10)), nrow = 5, ncol = 4,
                     dimnames = list(NULL, c("gen_1", "gen_2", "gen_3", "gen_4")))
    # 5 * 2*0.5*(1-0.5) = 2.5 for first column
    # 5 * 2*0.1*(1-0.1) = 5 * 0.18 = 0.9 for third/fourth
    result <- gt$sum_of_het_by_t(pmat)
    expect_equal(result$sum_of_het[1], 2.5)
    expect_equal(result$sum_of_het[2], 2.5)
    expect_equal(result$sum_of_het[3], 0.9)
    expect_equal(result$sum_of_het[4], 0.9)
  })
})

describe("conv_cor_wn_env", {
  it("returns a length-2 vector: numerator and denominator", {
    set.seed(42)
    pdiff <- matrix(rnorm(30), nrow = 5, ncol = 6)
    result <- gt$conv_cor_wn_env(pdiff)
    expect_equal(length(result), 2)
  })

  it("returns numerator smaller than denominator", {
    set.seed(99)
    pdiff <- matrix(rnorm(60), nrow = 10, ncol = 6)
    result <- gt$conv_cor_wn_env(pdiff)
    expect_true(abs(result[1]) <= result[2])
  })

  it("works with matrices of varying sizes", {
    set.seed(7)
    pdiff_4x3 <- matrix(rnorm(12), nrow = 4, ncol = 3)
    result <- gt$conv_cor_wn_env(pdiff_4x3)
    expect_equal(length(result), 2)
  })

  it("returns positive numerator for positively correlated data", {
    set.seed(123)
    a <- rnorm(10, mean = 0, sd = 1)
    b <- a + rnorm(10, mean = 0, sd = 0.1)  # strong correlation
    pdiff <- cbind(a, b)
    result <- gt$conv_cor_wn_env(pdiff)
    expect_true(result[1] > 0)
  })
  
  it("returns negative numerator for negatively correlated data", {
    set.seed(123)
    a <- rnorm(10, mean = 0, sd = 1)
    b <- -1*a + rnorm(10, mean = 0, sd = 0.1)  # strong correlation
    pdiff <- cbind(a, b)
    result <- gt$conv_cor_wn_env(pdiff)
    expect_true(result[1] < 0)
  })

  it("requires matrix input with ncol > 1 and nrow > 1", {
    expect_error(gt$conv_cor_wn_env(c(1, 2, 3)))
  })
})
