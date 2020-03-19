
library(stringr)
context("Naming cutpoints")
set.seed(1108)
x <- sample(0:4, size = 1000, replace = TRUE)
c1 <- c(1, 4)
c2 <- c(1, 2, 4)
c3 <- c(1, 2, 3, 4)

test_that("name_cuts outputs a character string", {
  expect_is(name_cuts(x, c1), "character")
})

test_that("name_cuts shows correct number of cuts in name", {
  expect_equal(str_count(name_cuts(x, c1), "v"), 1)
  expect_equal(str_count(name_cuts(x, c2), "v"), 2)
  expect_equal(str_count(name_cuts(x, c3), "v"), 3)
})

test_that("name_cuts has correct prefix based on number of intervals", {
  expect_true(grepl("^b", name_cuts(x, c1)))
  expect_true(grepl("^t", name_cuts(x, c2)))
  expect_true(grepl("^qd", name_cuts(x, c3)))
})
