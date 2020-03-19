
context("Parse Rd files")

pr <- parse_rd()
prs <- parse_rd(path = "function_dictionary.csv")

test_that("parse_rd has expected number of functions", {
  expect_equal(nrow(pr), length(list.files("man")))
})

test_that("table can be saved to file", {
  expect_identical(pr, prs)
  file.remove("function_dictionary.csv")
})

test_that("test", {
  expect_error(pr, NA)
  expect_error(prs, NA)
})
