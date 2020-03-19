
context("Find cutpoint by coxph")

data(lung)

test_that("findCutpointByCoxph works for univariable case", {
  suppressWarnings(obj <- findCutpointByCoxph(lung, Surv(time, status) ~ age))
  expect_length(obj, 5)
})

test_that("error if formula is multivariable", {
  expect_error(findCutpointByCoxph(lung, Surv(time, status) ~ age + sex))
})

test_that("error if explanatory variable is not continuous", {
  lung$wt.loss <- as.character(lung$wt.loss)
  expect_error(findCutpointByCoxph(lung, Surv(time, status) ~ wt.loss))
})

test_that("error if no variation in explanatory variable", {
  lung$dummy <- 1
  expect_error(findCutpointByCoxph(lung, Surv(time, status) ~ dummy))
})

test_that("unexpected error is thrown when there's numerical difficulty", {
  lung$time2 <- Inf
  expect_output(findCutpointByCoxph(lung, Surv(time2, status) ~ sex))
})
