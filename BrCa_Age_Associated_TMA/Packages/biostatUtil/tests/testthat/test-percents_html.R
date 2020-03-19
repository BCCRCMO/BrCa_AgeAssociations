
context("HTML output of Row and Column percentages")

library(htmlTable)
set.seed(13)
B <- matrix(rbinom(20, size = 20, prob = 0.3), nrow = 5,
            dimnames = list(paste0("Row", 1:5), paste0("Col", 1:4)))

test_that("HTML output works", {
  expect_error(htmlTable(colPercentAsHTML(
    B, keep = TRUE, caption = "Test Caption")),
    NA
    )
  expect_error(htmlTable(colPercentAsHTML(
    B, keep = FALSE, caption = "Test Caption")),
    NA
  )
  expect_error(htmlTable(rowPercentAsHTML(
    B, keep = TRUE, caption = "Test Caption")),
    NA
  )
  expect_error(htmlTable(rowPercentAsHTML(
    B, keep = FALSE, caption = "Test Caption")),
    NA
  )
  expect_error(htmlTable(rowColPercentAsHTML(
    B, keep = TRUE, caption = "Test Caption")),
    NA
  )
  expect_error(htmlTable(rowColPercentAsHTML(
    B, keep = FALSE, caption = "Test Caption")),
    NA
  )
})

test_that("Banded rows alternate per row group", {
  expect_error(htmlTable(colPercentAsHTML(
    B, keep = TRUE, banded.rows = TRUE, caption = "Test Caption")),
    NA
  )
  expect_error(htmlTable(colPercentAsHTML(
    B, keep = FALSE, banded.rows = TRUE, caption = "Test Caption")),
    NA
  )
  expect_error(htmlTable(rowPercentAsHTML(
    B, keep = TRUE, banded.rows = TRUE, caption = "Test Caption")),
    NA
  )
  expect_error(htmlTable(rowPercentAsHTML(
    B, keep = FALSE, banded.rows = TRUE, caption = "Test Caption")),
    NA
  )
  expect_error(htmlTable(rowColPercentAsHTML(
    B, keep = TRUE, banded.rows = TRUE, caption = "Test Caption")),
    NA
  )
  expect_error(htmlTable(rowColPercentAsHTML(
    B, keep = FALSE, banded.rows = TRUE, caption = "Test Caption")),
    NA
  )
})

test_that("Other formatting options work", {
  expect_error(htmlTable(rowPercentAsHTML(B, caption = "Test Caption")), NA)
  expect_error(htmlTable(colPercentAsHTML(B, html.table.border = 10)), NA)
  expect_error(htmlTable(rowColPercentAsHTML(
    B,
    keep = FALSE,
    banded.rows = TRUE,
    col.odd = "yellow",
    col.even = "green"
  )),
  NA)
})

test_that("Row and column names can be overwritten", {
  expect_error(htmlTable(colPercentAsHTML(
    B, keep = TRUE, row.names = letters[1:5], column.names = LETTERS[1:4])),
    NA
  )
  expect_error(htmlTable(rowPercentAsHTML(
    B, keep = TRUE, row.names = letters[1:5], column.names = LETTERS[1:4])),
    NA
  )
  expect_error(htmlTable(rowColPercentAsHTML(
    B, keep = TRUE, row.names = letters[1:5], column.names = LETTERS[1:4])),
    NA
  )
})

test_that("Transpose can be useful", {
  set.seed(13)
  C <- matrix(rbinom(10, 10, 0.4), nrow = 1)
  D <- matrix(rbinom(10, 10, 0.4), ncol = 1)
  expect_error(htmlTable(rowPercentAsHTML(C, transpose = FALSE)), NA)
  expect_error(htmlTable(rowPercentAsHTML(C, transpose = TRUE)), NA)
  expect_error(htmlTable(colPercentAsHTML(D, transpose = TRUE)), NA)
})
