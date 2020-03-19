
context("String utility functions")

test_that("Square brackets and parentheses are escaped", {
  expect_identical("\\[index\\]", escapeForGrep("[index]"))
  expect_identical("\\(parentheses\\)", escapeForGrep("(parentheses)"))
})

test_that("Indices change if case is considered", {
  expect_equal(5, indexOf("Animals", "a", ignore.case = FALSE))
  expect_equal(c(1, 5), indexOf("Animals", "a", ignore.case = TRUE))
})

test_that("NA returned if pattern longer than string or not found", {
  expect_equal(NA, indexOf("Animals", "Animalss"))
  expect_equal(NA, indexOf("Animals", "o"))
})

test_that("first letter lowercase is FALSE", {
  expect_false(isFirstLetterUpperCase("fALSE"))
})

test_that("first letter uppercase is TRUE", {
  expect_true(isFirstLetterUpperCase("True"))
})

test_that("first character empty string returns TRUE", {
  expect_true(isFirstLetterUpperCase(""))
})

test_that("First or all words capitalized", {
  expect_identical(
    "Ovarian Cancer",
    simpleCap("ovarian cancer", first.only = FALSE)
  )
  expect_identical(
    "Ovarian cancer",
    simpleCap("ovarian cancer", first.only = TRUE)
  )
})
