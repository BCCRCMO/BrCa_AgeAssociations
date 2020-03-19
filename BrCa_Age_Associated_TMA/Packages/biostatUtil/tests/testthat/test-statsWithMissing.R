
context("Statistics with Missing Values")

x <- c(10:1)
y <- c(1:10)
z <- c(10:1)
m <- c(x[1:9], NA)

test_that("Computations with missing codes work", {
  expect_equal(min(z[!(z %in% c(1, 3))]), minWithMissing(z, c(1, 3)))
  expect_equal(max(z[!(z %in% c(9, 10))]), maxWithMissing(z, c(9, 10)))
  expect_equal(mean(c(1:2, 8:10)), meanWithMissing(z, 3:7))
  expect_equal(sum(1:8), sumWithMissing(z, c(9, 10)))
  expect_equal(prod(1:8), prodWithMissing(z, c(9, 10)))

  r <- x / y
  r[(x %in% 1:2 | y %in% 1:2)] <- -1
  expect_equal(r, ratioWithMissing(x, y, 1:2))
})

test_that("Computations return missing code if all data is missing", {
  expect_equal(-1, minWithMissing(z, c(1:10)))
  expect_equal(-1, maxWithMissing(z, c(1:10)))
  expect_equal(-1, meanWithMissing(z, c(1:10)))
  expect_equal(-1, sumWithMissing(z, c(1:10)))
  expect_equal(-1, prodWithMissing(z, c(1:10)))
})
