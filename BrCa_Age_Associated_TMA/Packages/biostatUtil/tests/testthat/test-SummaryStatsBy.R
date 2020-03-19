
context("Generating cohort characteristics")
mtcars$vs <- as.factor(mtcars$vs)
mtcars$am <- as.factor(mtcars$am)
mtcars$onelev <- factor("d", levels = c("d", "e"))

test_that("No output if requested statistic not found", {
  stats.in1 <- c("median", "sd", "range", "mean", "missing", "cor")
  stats.in2 <- c("median", "sd", "range", "mean", "missing")
  ssb1 <-
    SummaryStatsBy(
      mtcars,
      by1 = "cyl",
      by2 = "gear",
      var.names = "mpg",
      stats = stats.in1
    )
  ssb2 <-
    SummaryStatsBy(
      mtcars,
      by1 = "cyl",
      by2 = "gear",
      var.names = "mpg",
      stats = stats.in2
    )
  expect_identical(ssb1, ssb2)
})

test_that("Order of input/output statistics match", {
  stats.in <- c("median", "range", "mean", "missing", "IQR")
  ssb <-
    SummaryStatsBy(
      mtcars,
      by1 = "cyl",
      by2 = "gear",
      var.names = "mpg",
      stats = stats.in
    )
  stats.out <- grep("\\*", rownames(ssb), value = TRUE, invert = TRUE)
  expect_identical(stats.in, stats.out)
})

test_that("Order of input/output statistics match if sd included", {
  stats.in <- c("median", "sd", "range", "mean", "missing", "IQR")
  ssb <-
    SummaryStatsBy(
      mtcars,
      by1 = "cyl",
      by2 = "gear",
      var.names = "mpg",
      stats = stats.in
    )
  stats.out <- grep("\\*", rownames(ssb), value = TRUE, invert = TRUE)
  expect_identical(grep("sd", stats.in, value = TRUE, invert = TRUE), stats.out)
})

test_that("Order of input/output variables match", {
  vars.in <- c("hp", "carb", "wt", "vs", "am", "qsec", "mpg")
  ssb <- SummaryStatsBy(mtcars, by1 = "cyl", by2 = "gear", var.names = vars.in)
  vars.out <- gsub("\\*", "", grep("\\*", rownames(ssb), value = TRUE))
  expect_identical(vars.in, vars.out)
})

test_that("Error if input isn't num, int, fac, char", {
  mtcars$logical_var <- TRUE
  expect_error(SummaryStatsBy(mtcars, by1 = "cyl", by2 = "gear",
                              var.names = "logical_var"))
})

test_that("Error thrown if argument names misspelled", {
  expect_error(
    SummaryStatsBy(mtcars, by1 = "cyl", by2 = "gear", var.names = c("mpgg"))
  )
  expect_error(
    SummaryStatsBy(mtcars, by1 = "cl", by2 = "ger", var.names = c("mpg"))
  )
})

test_that("Error if trying to split by one level factor", {
  expect_error(
    SummaryStatsBy(mtcars, by1 = "onelev", by2 = "gear", var.names = "mpg")
  )
  expect_error(
    SummaryStatsBy(mtcars, by1 = "gear", by2 = "onelev", var.names = "mpg")
  )
})

test_that("No error if factor has only one level", {
  expect_error(
    SummaryStatsBy(mtcars, by1 = "cyl", by2 = "gear", var.names = "onelev"),
    NA
  )
})

test_that("No errors for single continuous/categorical variable", {
  expect_error(
    SummaryStatsBy(mtcars, by1 = "cyl", by2 = "gear", var.names = c("mpg")),
    NA
  )
  expect_error(
    SummaryStatsBy(mtcars, by1 = "cyl", by2 = "gear", var.names = c("vs")),
    NA
  )
})

test_that("All formats work", {
  expect_error(
    SummaryStatsBy(mtcars, by1 = "cyl", by2 = "gear", var.names = "mpg", format = "raw"),
    NA
  )
  expect_error(
    SummaryStatsBy(mtcars, by1 = "cyl", by2 = "gear", var.names = "mpg", format = "pandoc"),
    NA
  )
  expect_error(
    SummaryStatsBy(mtcars, by1 = "cyl", by2 = "gear", var.names = "mpg", format = "html"),
    NA
  )
  expect_error(
    SummaryStatsBy(mtcars, by1 = "cyl", by2 = "gear", var.names = "mpg", format = "long"),
    NA
  )
})
