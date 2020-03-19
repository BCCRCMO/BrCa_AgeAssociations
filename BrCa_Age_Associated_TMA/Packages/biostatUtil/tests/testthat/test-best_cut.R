
context("Best cutpoint")

library(survival)
set.seed(1108)
n <- 4
lung_mod <- cbind(lung, fac = sample(1:n, nrow(lung), replace = TRUE))

test_that("factor with n levels can be separated into at most n groups", {
  expect_error(best_cut(
    Surv(time, status) ~ fac,
    lung_mod,
    n = "b",
    plot = FALSE
  ), NA)
  expect_error(best_cut(
    Surv(time, status) ~ fac,
    lung_mod,
    n = "t",
    plot = FALSE
  ), NA)
  expect_error(best_cut(
    Surv(time, status) ~ fac,
    lung_mod,
    n = "qd",
    plot = FALSE
  ), NA)
  expect_error(best_cut(
    Surv(time, status) ~ fac,
    lung_mod,
    n = "qn",
    plot = FALSE
  ))
})

test_that("survival curves can be saved or not saved", {
  expect_error(best_cut(
    Surv(time, status) ~ fac,
    lung_mod,
    n = "b",
    plot = TRUE
  ), NA)
  expect_error(best_cut(
    Surv(time, status) ~ fac,
    lung_mod,
    n = "b",
    plot = TRUE,
    filename = "test.png"
  ),
  NA)
  file.remove("test.png")
})

test_that("cutpoint criteria changes when likelihood is flat", {
  expect_error(best_cut(
    Surv(time, status) ~ fac,
    lung_mod,
    n = "b",
    plot = FALSE,
    AIC.range = 1e-2
  ),
  NA)
})

test_that("best cutpoint has lowest AIC when likelihood is not flat", {
  bc <- best_cut(Surv(time, status) ~ fac, lung_mod, n = "b", plot = FALSE)
  expect_true(which.min(unlist(bc$results$AIC)) == bc$opt.cut)
})
