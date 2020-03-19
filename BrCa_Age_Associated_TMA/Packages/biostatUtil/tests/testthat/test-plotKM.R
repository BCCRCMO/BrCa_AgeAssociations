context("Plot KM")

library(survival)
data(lung)

test_that("doKMPlots calls ggkm or plotKM", {
  expect_error(
    doKMPlots(lung, "time", "status", "sex", "Sex", use.ggkm = FALSE),
    NA
  )
  expect_error(
    doKMPlots(lung,
              "time",
              "status",
              "sex",
              "Sex",
              use.ggkm = TRUE,
              timeby = 200),
    NA
  )
})

test_that("factor variable drops unused levels", {
  lung$sex <- as.factor(lung$sex)
  expect_error(doKMPlots(lung, "time", "status", "sex", "Sex"), NA)
})

test_that("unused factor levels are removed", {
  lung$test <- factor(1:2, levels = 1:3)
  expect_error(
    doKMPlots(lung, "time", "status", "test", "Test", use.ggkm = FALSE),
    NA
  )
})

test_that("other plotting options can be specified", {
  expect_error(doKMPlots(lung, "time", "status", "sex", "Sex", use.ggkm = FALSE,
                         grey.scale = TRUE, shading.colors = NULL,
                         legend.pos = "top"), NA)
  expect_error(doKMPlots(lung, "time", "status", "sex", "Sex", use.ggkm = TRUE,
                         timeby = 200, cox.ref.group = "2", show.risk = FALSE,
                         use.firth = 0.8), NA)
})

test_that("plot can be saved to file", {
  expect_error(doKMPlots(lung, "time", "status", "sex", "Sex", use.ggkm = FALSE,
                         file.name = "test.pdf"), NA)
  expect_error(doKMPlots(lung, "time", "status", "sex", "Sex", use.ggkm = FALSE,
                         file.name = "test.odf"))
  file.remove(list.files(pattern = "test\\."))
})

test_that("plot statistics indicated by reference group", {
  expect_error(doKMPlots(lung, "time", "status", "sex", "Sex", use.ggkm = FALSE,
                         km.plot.ref.group = "2"), NA)
})

test_that("survival fits can be compared", {
  p1 <- doKMPlots(lung, "time", "status", "sex", "Sex", use.ggkm = TRUE,
                  timeby = 200, use.firth = 0.8, cox.ref.group = "2")
  p2 <- doKMPlots(lung, "time", "status", "sex", "Sex", use.ggkm = TRUE,
                  timeby = 200, use.firth = 0.8, cox.ref.group = "2",
                  sfit2 = survfit(Surv(time, status) ~ age, lung))
  p3 <- doKMPlots(lung, "time", "status", "sex", "Sex", use.ggkm = TRUE,
                  timeby = 200, use.firth = -1, cox.ref.group = "2",
                  sfit2 = survfit(Surv(time, status) ~ age, lung))
  expect_false(identical(p1, p2))
  expect_false(identical(p2, p3))
})

test_that("margin widths adapt to length of predictor variable labels", {
  lung$gender <-
    factor(lung$sex,
           levels = c(1, 2),
           labels = c("Male", "Female"))
  lung$gender2 <-
    factor(lung$sex,
           levels = c(1, 2),
           labels = c("Male", "Femal"))
  p4 <- doKMPlots(lung, "time", "status", "gender", "Sex", use.ggkm = TRUE,
                  timeby = 200, use.firth = 0.8, cox.ref.group = "2")
  p5 <- doKMPlots(lung, "time", "status", "gender2", "Sex", use.ggkm = TRUE,
                  timeby = 200, use.firth = 0.8, cox.ref.group = "2")
  expect_error(p4, NA)
  expect_error(p5, NA)
})

test_that("HR is shown for all levels of multilevel predictor", {
  lung$ph.ecog.f <- factor(lung$ph.ecog)
  p6 <- doKMPlots(lung, "time", "status", "ph.ecog.f", "PH.ECOG",
                  use.ggkm = TRUE, timeby = 200, use.firth = 0.8,
                  cox.ref.group = "0")
  expect_error(p6, NA)
})
