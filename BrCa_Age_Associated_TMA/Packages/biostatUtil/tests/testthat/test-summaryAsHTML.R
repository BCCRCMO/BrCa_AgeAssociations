
context("summary as HTML format")

library(htmlTable)
set.seed(1)
x <- rnorm(100)

test_that("Numerical summary as HTML format works", {
  expect_error(htmlTable(summaryAsHTML(x)), NA)
})
