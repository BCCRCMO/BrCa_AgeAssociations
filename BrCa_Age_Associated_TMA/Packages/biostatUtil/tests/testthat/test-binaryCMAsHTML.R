
context("Binary class confusion matrices as HTML format")

library(htmlTable)
set.seed(547)
n <- 20
x <- rbinom(n, size = 1, prob = 0.6)
y <- rbinom(n, size = 1, prob = 0.4)
results <- binaryCMAsHTML(x, y, "Test", num.boot = 1000)

test_that("binaryCMAsHTML works", {
  expect_error(htmlTable(results), NA)
})
