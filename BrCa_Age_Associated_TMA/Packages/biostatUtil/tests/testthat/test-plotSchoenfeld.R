
context("Plot Schoenfeld residual plot")

library(survival)
test1 <- list(time = c(4, 3, 1, 1, 2, 2, 3),
              status = c(1, 1, 1, 0, 1, 1, 0),
              x = c(0, 2, 1, 1, 1, 0, 0),
              sex = c(0, 0, 0, 0, 1, 1, 1))

test_that("Schoenfeld residual plot works whether variable to plot is given or not", {
  expect_error(
    plotSchoenfeld(test1, Surv(time, status) ~ x + strata(sex), "x"),
    NA
  )
  expect_error(
    plotSchoenfeld(test1, Surv(time, status) ~ x + strata(sex)),
    NA
  )
})
