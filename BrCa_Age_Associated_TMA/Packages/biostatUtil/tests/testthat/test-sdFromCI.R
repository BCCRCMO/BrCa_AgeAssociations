
context("Standard deviation from Confidence Interval")

test_that("function returns a list of length two named 'sd' and 'Pval'", {
  x <- sdFromCI(0.43, 0.13, 1.37)
  expect_is(x, "list")
  expect_length(x, 2)
  expect_identical(names(x), c("sd", "Pval"))
})
