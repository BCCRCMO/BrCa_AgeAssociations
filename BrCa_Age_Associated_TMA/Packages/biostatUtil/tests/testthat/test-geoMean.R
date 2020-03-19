
context("Geometric mean")

set.seed(1)
x <- rexp(100, rate = 0.5)
y <- rexp(100, rate = 0.3)
y[3] <- NA
z <- rep(NA, 100)

test_that("Geometric mean works", {
  expect_is(geoMean(x), "numeric")
  expect_is(geoMean(x, na.rm = TRUE), "numeric")
})

test_that("Missing elements are handled", {
  expect_is(geoMean(y, na.rm = TRUE), "numeric")
  expect_equal(NA_real_, geoMean(y))
  expect_equal(NA, geoMean(z, na.rm = TRUE))
})
