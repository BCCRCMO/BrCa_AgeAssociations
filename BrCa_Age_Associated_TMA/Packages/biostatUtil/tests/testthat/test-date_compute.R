
context("Date computation functions")

x5 <- "1992/01/27"
x6 <- "2003/03/21"

test_that("Adding date works for allowed units", {
  expect_identical(
    "1992/02/06",
    addToDate(x5, 10, date.format = "YYYY.MM.DD", units = "days")
  )
  expect_identical(
    "1992/04/06",
    addToDate(x5, 10, date.format = "YYYY.MM.DD", units = "weeks")
  )
  expect_identical(
    "1992/11/27",
    addToDate(x5, 10, date.format = "YYYY.MM.DD", units = "months")
  )
  expect_identical(
    "2002/01/26",
    addToDate(x5, 10, date.format = "YYYY.MM.DD", units = "years")
  )
})

test_that("Adding date handles erroneous inputs", {
  expect_equal(NA, addToDate(NA, 10))
  expect_equal(NA, addToDate(x5, NA))
  expect_equal(
    "Missing",
    addToDate(
      x5,
      10,
      existing.missing.codes = x5,
      return.missing.code = "Missing"
    )
  )
})

test_that("Subtracting date works for allowed units", {
  expect_equal(
    4071,
    diffDate(x6, x5, date.format = "YYYY.MM.DD", units = "days")
  )
  expect_equal(
    581.5714,
    diffDate(x6, x5, date.format = "YYYY.MM.DD", units = "weeks"),
    tolerance = .002
  )
  expect_equal(
    133.4754,
    diffDate(x6, x5, date.format = "YYYY.MM.DD", units = "months"),
    tolerance = .002
  )
  expect_equal(
    11.14607,
    diffDate(x6, x5, date.format = "YYYY.MM.DD", units = "years"),
    tolerance = .002
  )
})

test_that("Subtracting date handles erroneous inputs", {
  expect_equal(NA, diffDate(x6, NA))
  expect_equal(NA, diffDate(NA, x5))
  expect_equal(
    "Missing",
    diffDate(
      x6,
      x5,
      existing.missing.codes = x5,
      return.missing.code = "Missing"
    )
  )
  expect_equal(
    "Missing",
    diffDate(
      x6,
      x5,
      existing.missing.codes = x6,
      return.missing.code = "Missing"
    )
  )
})

test_that("Comparing dates works for any two valid dates", {
  expect_equal(
    -1,
    compareDate("01/22/1949", "04/13/1950", date.format = "MM.DD.YYYY")
  )
  expect_equal(
    0,
    compareDate("04/13/1950", "04/13/1950", date.format = "MM.DD.YYYY")
  )
  expect_equal(
    1,
    compareDate("04/13/1959", "04/13/1950", date.format = "MM.DD.YYYY")
  )
  expect_equal(
    NA,
    compareDate(NA, "04/13/1950", date.format = "MM.DD.YYYY")
  )
})
