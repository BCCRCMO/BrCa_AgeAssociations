
context("Cox proportional hazards models utility functions")

library(coxphf)
test1 <- list(time = c(4, 3, 1, 1, 2, 2, 3),
              status = c(1, 1, 1, 0, 1, 1, 0),
              x = c(0, 2, 1, 1, 1, 0, 0),
              sex = c(0, 0, 0, 0, 1, 1, 1))
mod1 <- coxph(Surv(time, status) ~ x + strata(sex), test1)

test2 <- data.frame(list(
  start = c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
  stop = c(2, 3, 6, 7, 8, 9, 9, 9, 14, 17),
  event = c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
  x    = c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0)
))
mod2 <- coxphf(formula = Surv(start, stop, event) ~ x, pl = FALSE, data = test2)

test_that("coxphOut matrix summary, with overwritable rownames", {
  expect_is(coxphOut(mod1), "matrix")
  expect_length(coxphOut(mod1), 9)
  expect_identical("x", rownames(coxphOut(mod1)))
  expect_identical("y", rownames(coxphOut(mod1, coefnames = "y")))
})

test_that("Xunivcoxph prints HR & CI for coxph/coxphf, but nothing else", {
  x1 <- coxphOut(mod1, conf.level = 0.95, digits = 2)
  x1 <- paste0("HR ", x1[, "HR"], " (95% CI: ",
               x1[, "2.5 %"], "-", x1[, "97.5 %"], ")")
  y1 <- Xunivcoxph(mod1, digits = 2)
  expect_identical(x1, y1)

  x2 <- coxphOut(mod2, conf.level = 0.95, digits = 2)
  x2 <- paste0("HR(F) ", x2[, "HR"], " (95% CI: ",
               x2[, "2.5 %"], "-", x2[, "97.5 %"], ")")
  y2 <- Xunivcoxph(mod2, digits = 2)
  expect_identical(x2, y2)

  expect_error(
    Xunivcoxph(survdiff(Surv(time, status) ~ x + strata(sex), test1))
  )
})

test_that("printCoxMod prints HTML Cox model output", {
  expect_output(printCoxMod(summary(mod1)$coef))
})
