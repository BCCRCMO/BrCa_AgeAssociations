
context("Round small numbers")

test_that("vector input works", {
  expect_error(round_small(2e-04), NA)
})

test_that("matrix input works", {
  set.seed(12)
  x <- matrix(rexp(25, 3), nrow = 5)
  expect_error(round_small(x, digits = 1), NA)
})

test_that("scientific notation works", {
  expect_true(grepl("e", round_small(2e-4, sci = TRUE)))
})

test_that("NA doesn't affect rounding", {
  expect_error(round_small(c(2e-05, 3e-04, NA, 4e-02)), NA)
})
