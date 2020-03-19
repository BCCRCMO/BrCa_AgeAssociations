
context("Bootstrapped mean")

mySeed <- 789
set.seed(mySeed)
test.vec <- rnorm(n = 60, mean = 7, sd = 3)
test.mat <- matrix(test.vec, nrow = 10)
vecTestResult <- bootMean(test.vec)
matTestResult <- bootMean(test.mat)

set.seed(mySeed)
s <- ifelse(rexp(100, 0.5) < 1, NA, rexp(100, 0.5))
sTestResultNaRm <- bootMean(s, na.rm = TRUE)

test_that("result should be equal regardless of input class", {
  expect_identical(vecTestResult, matTestResult)
})

test_that("vector output has correct values based on random seed", {
  expect_equal(vecTestResult$n, length(test.vec))
  expect_equal(length(vecTestResult$obs.mean), 1)
  expect_equal(signif(vecTestResult$obs.mean), signif(6.725579))
  expect_equal(signif(vecTestResult$ci[1]), signif(5.983460))
  expect_equal(signif(vecTestResult$ci[2]), signif(7.381294))
})

test_that("matrix output has correct values based on random seed", {
  expect_equal(matTestResult$n, length(test.mat))
  expect_equal(length(matTestResult$obs.mean), 1)
  expect_equal(signif(matTestResult$obs.mean), signif(6.725579))
  expect_equal(signif(matTestResult$ci[1]), signif(5.983460))
  expect_equal(signif(matTestResult$ci[2]), signif(7.381294))
})

test_that("NAs are handled when na.rm = TRUE", {
  expect_equal(signif(sTestResultNaRm$obs.mean), signif(2.195333))
  expect_equal(signif(sTestResultNaRm$n), length(s))
  expect_equal(signif(sTestResultNaRm$ci[1]), signif(1.659419))
  expect_equal(signif(sTestResultNaRm$ci[2]), signif(2.715013))
})

test_that("warnings if input is a character", {
  expect_warning(bootMean("xx"))
})

test_that("NAs not handled when na.rm = FALSE", {
  expect_equal(bootMean(s)$obs.mean, NA_real_)
})
