
context("Explore Data")

test_that("PDF output generated", {
  mtcars$vs <- as.factor(mtcars$vs)
  mtcars$am <- as.factor(mtcars$am)
  exploreData(mtcars)
  expect_true("DataSummary.pdf" %in% list.files())
  file.remove("DataSummary.pdf")
})
