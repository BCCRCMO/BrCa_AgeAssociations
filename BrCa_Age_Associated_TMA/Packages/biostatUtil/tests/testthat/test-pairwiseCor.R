
context("Pairwise correlations")

set.seed(123)
x <- data.frame(matrix(rnorm(25), nrow = 5))
y <- pairwiseCor(x)

test_that("Number of rows equals number of distinct pairs", {
  expect_equal(nrow(y), choose(ncol(x), 2))
})

test_that("Doesn't work when some columns aren't numeric", {
  x[, 4] <- letters[1:5]
  expect_error(pairwiseCor(x))
})
