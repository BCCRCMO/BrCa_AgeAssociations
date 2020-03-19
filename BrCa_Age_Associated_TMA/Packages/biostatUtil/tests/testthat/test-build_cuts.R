
context("Building cutpoints")
set.seed(1108)
x <- sample(0:4, size = 1000, replace = TRUE)

test_that("build_cuts outputs correct data structure", {
  expect_is(build_cuts(x, list = TRUE), "list")
  expect_is(build_cuts(x), "data.frame")
})

test_that(
  "build_cuts gives choose(dplyr::n_distinct(x) - 1, n - 1) elements/columns", {
    expect_equal(ncol(build_cuts(x, n = "t")), choose(4, 2))
    expect_equal(ncol(build_cuts(x, n = "qd")), choose(4, 3))
    expect_length(build_cuts(x, n = "qd", list = TRUE), choose(4, 3))
  }
)

test_that("prefix added in variable name", {
  expect_true(all(grepl("gene1_", names(build_cuts(x, var.prefix = "gene1")))))
})

test_that("error if cutting character vector", {
  expect_error(build_cuts(sample(LETTERS[1:3], size = 100, replace = TRUE)))
})

test_that("error if cutting into more than one number of intervals", {
  expect_error(build_cuts(x, n = c("b", "t")))
})

test_that("error if number of intervals exceeds number of levels", {
  expect_error(build_cuts(sample(0:3, size = 1000, replace = TRUE), n = "qn"))
})
