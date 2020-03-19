
context("G-Rho Tests for Survival Differences")
library(survival)

test_that("formula and data work", {
  expect_error(grhoTests(Surv(futime, fustat) ~ rx, ovarian, digits = 2), NA)
})
