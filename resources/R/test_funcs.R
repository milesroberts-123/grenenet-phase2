library(testthat)

# Source the functions under test
box::use(./misc)
box::use(./gt)
box::use(./plot)

describe("fitfreq", {
  it("returns 0.5 when q=0.5, h=0.5, s=0", {
    expect_equal(misc$fitfreq(0.5, 0.5, 0), 0.5)
  })

  it("returns the same allele frequency when s=0 regardless of h", {
    expect_equal(misc$fitfreq(0.2, 0.0, 0), 0.2)
    expect_equal(misc$fitfreq(0.8, 1.0, 0), 0.8)
  })

  it("returns q when h=0 and s!=0 (dominance coefficient is zero)", {
    expect_equal(misc$fitfreq(0.3, 0, 0.1), 0.3)
  })

  it("throws error or returns valid value for q in [0,1]", {
    expect_equal(misc$fitfreq(0.0, 0.5, 0.1), 0.0)
    expect_equal(misc$fitfreq(1.0, 0.5, 0.1), (1 + 0.1) / (1 + 0.1))
  })
})

describe("freq_increments", {
  it("correctly computes differences between adjacent generations", {
    pmat <- matrix(c(0.1, 0.2, 0.3, 0.4,
                     0.2, 0.3, 0.4, 0.5), nrow = 2, ncol = 4)
    expected <- matrix(c(0.1, 0.1, 0.1,
                         0.1, 0.1, 0.1), nrow = 2, ncol = 3)
    result <- gt$freq_increments(pmat)
    expect_true(all(result == expected))
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
    expect_equal(result, matrix(c(0.2, 0.3), nrow = 1))
  })
})

describe("ssh", {
  it("correctly computes sum of 2*x*(1-x) for a vector", {
    expect_equal(gt$ssh(c(0.5)), 0.5)
    expect_equal(gt$ssh(c(0.0, 1.0)), 0)
    expect_equal(gt$ssh(c(1/3, 1/3, 1/3)), 2 * (1/3) * (2/3) * 3)
  })

  it("returns zero for fixation states (all 0 or all 1)", {
    expect_equal(gt$ssh(c(0, 0, 0)), 0)
    expect_equal(gt$ssh(c(1, 1, 1)), 0)
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
    expect_true(result[1, 3] == 6)
    expect_true(result[2, 3] == 9)
    expect_true(result[1, 2] == 2)
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
    covmat <- matrix(c(4, 1, 0.5, 9), nrow = 2)
    result <- gt$gt_from_covmat(covmat)
    expect_equal(result$sum_var[2], 13)  # 4+9
    expect_equal(result$sum_cov[2], 1)
    expect_equal(result$sum_abs_cov[2], 2)
  })
})
