
context("Compute cohort characteristics")

library(htmlTable)
mtcars$vs <- as.factor(mtcars$vs)

test_that("defaults can be changed and output is a length 4 list", {
  dcc <- suppressWarnings(doCohortCharacteristics(
    input.d = mtcars,
    marker.name = "vs",
    marker.description = "V/S",
    var.names = c("disp", "hp", "qsec"),
    var.descriptions = c("displacement", "horsepower", "1/4 mile time"),
    is.var.continuous = c(TRUE, TRUE, TRUE),
    do.droplevels = TRUE,
    caption = "Some mtcars summaries",
    custom.total.label = "TOTAL",
    custom.marker.labels = c("V/S = 0", "V/S = 1"),
    stat.tests = c("spearman", "kruskal", "wilcox")))
  expect_length(dcc, 4)
  expect_is(dcc$result.table, "matrix")
  expect_length(dcc$stat.tests.results, 3)
  expect_is(dcc$result.table.html, "character")
  expect_is(dcc$result.table.bamboo, "character")
})

test_that("statistical tests for categorical variables work", {
  dcc_kendall <- suppressWarnings(doCohortCharacteristics(
    input.d = mtcars, marker.name = "vs", marker.description = "cylinders",
    var.names = "am", var.descriptions = c("transmission"),
    is.var.continuous = FALSE, caption = "Some mtcars summaries",
    stat.tests = "kendall"))
  dcc_chisq <- doCohortCharacteristics(
    input.d = mtcars, marker.name = "vs", marker.description = "cylinders",
    var.names = "am", var.descriptions = c("transmission"),
    is.var.continuous = FALSE, caption = "Some mtcars summaries",
    stat.tests = "chisq")
  dcc_fisher <- doCohortCharacteristics(
    input.d = mtcars, marker.name = "vs", marker.description = "cylinders",
    var.names = "am", var.descriptions = c("transmission"),
    is.var.continuous = FALSE, caption = "Some mtcars summaries",
    stat.tests = "fisher")

  mtcars$vs <- as.factor(mtcars$vs)
  mtcars$am <- as.factor(mtcars$am)
  dcc_cm_marker <- doCohortCharacteristics(
    input.d = mtcars, marker.name = "vs", marker.description = "cylinders",
    var.names = c("am"), var.descriptions = c("transmission"),
    is.var.continuous = c(FALSE), caption = "Some mtcars summaries",
    stat.tests = "confusionMarkerAsRef")
  dcc_cm_var <- doCohortCharacteristics(
    input.d = mtcars, marker.name = "vs", marker.description = "cylinders",
    var.names = c("am"), var.descriptions = c("transmission"),
    is.var.continuous = c(FALSE), caption = "Some mtcars summaries",
    stat.tests = "confusionVarAsRef")
  expect_length(dcc_kendall$stat.tests.results, 1)
  expect_length(dcc_chisq$stat.tests.results, 1)
  expect_length(dcc_fisher$stat.tests.results, 1)
  expect_error(dcc_cm_marker$stat.tests.results, NA)
  expect_error(dcc_cm_var$stat.tests.results, NA)
})

test_that("percentages can be for marginal rows or columns only", {
  dcc_row <- doCohortCharacteristics(
    input.d = mtcars,
    marker.name = "wt",
    marker.description = "weight",
    var.names = c("disp", "hp", "am"),
    var.descriptions = c("displacement", "horsepower", "transmission"),
    is.var.continuous = c(TRUE, TRUE, FALSE),
    show.percent = "row")
  dcc_col <- doCohortCharacteristics(
    input.d = mtcars,
    marker.name = "wt",
    marker.description = "weight",
    var.names = c("disp", "hp", "am"),
    var.descriptions = c("displacement", "horsepower", "transmission"),
    marker.value.labels.tolower = FALSE,
    is.var.continuous = c(TRUE, TRUE, FALSE),
    show.percent = "column")
  expect_error(dcc_row, NA)
  expect_error(dcc_col, NA)
})

test_that("error when marker not found", {
  expect_error(doCohortCharacteristics(mtcars, marker.name = "click",
                                       var.names = c("A", "B")))
})
