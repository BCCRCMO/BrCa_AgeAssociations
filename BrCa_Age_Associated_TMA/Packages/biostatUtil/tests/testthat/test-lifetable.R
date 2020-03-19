
context("Lifetables")

library(survival)
obj <- survfit(Surv(futime, fustat) ~ rx, data = ovarian)

test_that("it can print with or without summary", {
  expect_error(lifetable(obj), NA)
  expect_error(lifetable(obj, summary = TRUE), NA)
})

test_that("it works with one time interval", {
  expect_error(lifetable(obj, ntimes = 1), NA)
})

test_that("can hide strata name", {
  expect_error(lifetable(obj, show.strata = FALSE), NA)
})

test_that("dimensions are correct", {
  lto <- lifetable(obj, ntimes = 5)
  expect_length(lto, 12)
  expect_equal(nlevels(lto$strata) * nlevels(lto$times), nrow(lto))
})
