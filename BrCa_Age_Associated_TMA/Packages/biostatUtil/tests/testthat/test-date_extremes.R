
context("Date extremes")

test_that("Min and max work for vectors", {
  t1 <- as.Date(c("2015-03-01", "2015-02-15", NA, "2014-05-01"))
  expect_identical(as.Date("2015-03-01"), maxDate(t1))
  expect_identical(as.Date("2014-05-01"), minDate(t1))
})

test_that("Min and max work for arrays", {
  t2 <- c("03/21/1992", "04/21/2013", "10/10/2015")
  expect_identical("10/10/2015", maxDateArray(t2))
  expect_identical("03/21/1992", minDateArray(t2))
})

test_that("Array versions handle missing values", {
  expect_equal(NA, maxDateArray(character(0)))
  expect_equal(NA, minDateArray(character(0)))
  expect_identical(
    "Missing",
    maxDateArray(
      "09/11/2011",
      existing.missing.codes = "09/11/2011",
      return.missing.code = "Missing"
    )
  )
  expect_identical(
    "Missing",
    minDateArray(
      "09/11/2011",
      existing.missing.codes = "09/11/2011",
      return.missing.code = "Missing"
    )
  )
})
