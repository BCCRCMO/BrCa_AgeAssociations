
context("Date formatting functions")

x1 <- "2014-07-08"
x2 <- "2014/07/08"
x3 <- "2014|07|08"
x4 <- "20140708"

test_that(
  "Wrapper around numeric method for as.Date uses 1970 as default origin", {
    expect_identical(as.Date("1970-01-11"), numericToDate(10))
    expect_identical(as.Date("2000-09-21"), numericToDate(10, "2000-09-11"))
  }
)

test_that("Any date format works with any separator", {
  expect_identical("%m-%d-%Y", getFormat(x1, "MM.DD.YYYY"))
  expect_identical("%b-%d-%Y", getFormat(x1, "MMM.DD.YYYY"))
  expect_identical("%d-%m-%Y", getFormat(x1, "DD.MM.YYYY"))
  expect_identical("%d-%b-%Y", getFormat(x1, "DD.MMM.YYYY"))
  expect_identical("%Y-%m-%d", getFormat(x1, "YYYY.MM.DD"))
  expect_identical("%Y-%b-%d", getFormat(x1, "YYYY.MMM.DD"))

  expect_identical("%m/%d/%Y", getFormat(x2, "MM.DD.YYYY"))
  expect_identical("%b/%d/%Y", getFormat(x2, "MMM.DD.YYYY"))
  expect_identical("%d/%m/%Y", getFormat(x2, "DD.MM.YYYY"))
  expect_identical("%d/%b/%Y", getFormat(x2, "DD.MMM.YYYY"))
  expect_identical("%Y/%m/%d", getFormat(x2, "YYYY.MM.DD"))
  expect_identical("%Y/%b/%d", getFormat(x2, "YYYY.MMM.DD"))

  expect_identical("%m|%d|%Y", getFormat(x3, "MM.DD.YYYY"))
  expect_identical("%b|%d|%Y", getFormat(x3, "MMM.DD.YYYY"))
  expect_identical("%d|%m|%Y", getFormat(x3, "DD.MM.YYYY"))
  expect_identical("%d|%b|%Y", getFormat(x3, "DD.MMM.YYYY"))
  expect_identical("%Y|%m|%d", getFormat(x3, "YYYY.MM.DD"))
  expect_identical("%Y|%b|%d", getFormat(x3, "YYYY.MMM.DD"))
})

test_that("Custom separator can be added if '-', '/', '|' not used", {
  expect_identical("%m-%d-%Y", getFormat(x4, "MM.DD.YYYY", sep = "-"))
})

test_that("Error if invalid format provided", {
  expect_error(getFormat(x1, "MM.DD.FFFF"))
})

test_that("Format date works for any format", {
  expect_identical("07/08/2011", formatDate(8, 7, 2011))
  expect_identical(
    "Jan-10-2015",
    formatDate(10, 1, 2015, date.format = "MMM.DD.YYYY", sep = "-")
  )
  expect_identical(
    "08-07-2011",
    formatDate(8, 7, 2011, date.format = "DD.MM.YYYY", sep = "-")
  )
  expect_identical(
    "08|Jul|2011",
    formatDate(8, 7, 2011, date.format = "DD.MMM.YYYY", sep = "|")
  )
  expect_identical(
    "2011/07/08",
    formatDate(8, 7, 2011, date.format = "YYYY.MM.DD")
  )
  expect_identical(
    "2016/Sep/11",
    formatDate(11, 9, 2016, date.format = "YYYY.MMM.DD")
  )
})

test_that("Numeric dates can be reformatted for certain original formats", {
  expect_identical(
    "1991/Sep/11",
    cleanDate(
      11091991,
      original.format = "DD.MM.YYYY",
      preferred.format = "YYYY.MMM.DD"
    )
  )
  expect_identical(
    "1991/Sep/11",
    cleanDate(
      09111991,
      original.format = "MM.DD.YYYY",
      preferred.format = "YYYY.MMM.DD"
    )
  )
  expect_identical(
    "1991/Sep/11",
    cleanDate(
      19910911,
      original.format = "YYYY.MM.DD",
      preferred.format = "YYYY.MMM.DD"
    )
  )
  expect_error(cleanDate(19910911, original.format = "YYYY.MMM.DD"))
})

test_that("Character dates can be reformatted for certain original formats", {
  expect_identical(
    "1991/Sep/13",
    cleanDate(
      "13/09/1991",
      original.format = "DD.MM.YYYY",
      preferred.format = "YYYY.MMM.DD"
    )
  )
  expect_identical(
    "1991/Sep/13",
    cleanDate(
      "09/13/1991",
      original.format = "MM.DD.YYYY",
      preferred.format = "YYYY.MMM.DD"
    )
  )
  expect_identical(
    "1991/Sep/11",
    cleanDate(
      "09/11/1991",
      original.format = "MM.DD.YYYY",
      preferred.format = "YYYY.MMM.DD"
    )
  )
  expect_identical(
    "1991/Nov/09",
    cleanDate(
      "09/11/1991",
      original.format = "DD.MM.YYYY",
      preferred.format = "YYYY.MMM.DD"
    )
  )
  expect_identical(
    "1991/Sep/11",
    cleanDate(
      "1991/09/11",
      original.format = "YYYY.MM.DD",
      preferred.format = "YYYY.MMM.DD"
    )
  )
  expect_error(cleanDate("09/11/1991", original.format = "YYYY.MMM.DD"))
})

test_that("Reformatting dates handles erroneous inputs", {
  expect_equal("Missing", cleanDate(NA, return.missing.code = "Missing"))
  expect_equal(x4, cleanDate(x4, existing.missing.codes = x4))
  expect_equal(NA, cleanDate("This is not a date"))
})
