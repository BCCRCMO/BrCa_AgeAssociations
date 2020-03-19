
context("Table as HTML format")

library(htmlTable)
A <- matrix(c(2, 3, 5, 10), nrow = 2,
            dimnames = list(c("Row1", "Row2"), c("Col1", "Col2")))

test_that("row and colum names can be overridden", {
  expect_error(htmlTable(tableAsHTML(A, row.names = c("R1", "R2"),
                                     column.names = c("C1", "C2"))), NA)
})

test_that("matrix with banded rows works", {
  expect_error(htmlTable(tableAsHTML(A, banded.rows = FALSE)), NA)
  expect_error(htmlTable(tableAsHTML(A, banded.rows = TRUE)), NA)
})
