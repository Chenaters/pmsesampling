library(pmsesampling)
library(testthat)

test_that("pmse_samplesize returns a valid positive numeric sample size", {
  result = pmse_samplesize(k = 10, p = 5, sigma_k2 = 0.5, sigma_p2 = 0.6)


  expect_type(result, "double")
  expect_gt(result, 0)
})

test_that("pmse_samplesize returns a positive numeric sample size from sample covariance matrix", {
  file_path = file.path(testthat::test_path(), "testdata/SIGMA.txt")
  file_contents = Matrix::forceSymmetric(as.matrix(read.table(file_path)))
  SD = 1
  cov = sweep(sweep(file_contents, 1, SD, "*"), 2, SD, "*")
  result = pmse_samplesize(k = 12, p = 3, cov = cov)

  expect_type(result, "double")
  expect_equal(result, 103.6487, tolerance = 1e-6)

})

test_that("pmse_samplesize returns an error since it does not have predictor amount", {
  expect_error(pmse_samplesize(sigma_k2 = 0.5, sigma_p2 = 0.6))

})

test_that("pmse_samplesize returns an error since its unable to produce predictor error variance", {
  # For example, negative values might be invalid.
  expect_error(pmse_samplesize(k = 10, p = 3))
})

test_that("direct sigma inputs yield positive size", {
  out <- pmse_samplesize(k = 10, p = 5, sigma_k2 = 0.5, sigma_p2 = 0.6)
  expect_type(out, "double")
  expect_gt(out, 0)
})

test_that("covariance-matrix path matches known reference", {
  file_path = file.path(testthat::test_path(), "testdata/SIGMA.txt")
  SIG  <- Matrix::forceSymmetric(as.matrix(read.table(file_path)))
  size <- pmse_samplesize(k = 12, p = 3, cov = SIG)
  expect_equal(size, 103.6487, tolerance = 1e-4, scale = 100)
})

test_that("correlation-matrix path is accepted", {
  C  <- diag(1, 3); C[1,2] <- C[2,1] <- .3; SD <- c(1, 2, 3)
  n  <- pmse_samplesize(k = 3, p = 1, corr = C, SD = SD)
  expect_true(is.numeric(n) && n > 0)
})

test_that("f2 / R2 path returns efficient size only", {
  n <- pmse_samplesize(k = 4, p = 2, f2_2 = .15)
  expect_length(n, 1)
  expect_gt(n, 0)
})

test_that("function errors on unattainable PMSE", {
  expect_error(
    pmse_samplesize(k = 5, p = 2, sigma_k2 = 1, PMSE_val_k = 0.4)
  )
})
