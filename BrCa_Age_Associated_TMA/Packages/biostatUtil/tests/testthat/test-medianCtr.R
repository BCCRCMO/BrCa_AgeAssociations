
context("Median center")

test_that("median centering works", {
  set.seed(1)
  x <- matrix(rnorm(200), nrow = 10)
  expect_error(medianCtr(x), NA)
})

test_that("error thrown if input is NA", {
  expect_error(medianCtr(NA))
})
