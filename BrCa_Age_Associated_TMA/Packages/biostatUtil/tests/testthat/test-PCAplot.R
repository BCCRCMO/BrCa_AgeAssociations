
context("Principal component scatterplots and multiple plotting")

test_that("PCAplot works", {
  expect_error(PCAplot(mtcars$cyl, mtcars), NA)
})

test_that("One plot can still use multiplot", {
  p <- ggplot2::ggplot(mtcars, aes(wt, mpg)) +
    ggplot2::geom_point()
  expect_error(multiplot(p), NA)
})
