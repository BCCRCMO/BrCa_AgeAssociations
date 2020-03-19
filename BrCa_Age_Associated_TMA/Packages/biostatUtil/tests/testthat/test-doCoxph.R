context("Cox model generic")

library(survival)
data(lung)
set.seed(1)
lung$time2 <- sample(lung$time, replace = TRUE)
lung$status2 <- sample(lung$status, replace = TRUE)

test_that("doCoxphGeneric HRs are different when reference group changes", {
  res1 <- doCoxphGeneric(input.d = lung,
                         var.names = "sex",
                         var.descriptions = "Sex",
                         show.var.detail = TRUE,
                         var.names.surv.time = "time",
                         var.names.surv.status = "status",
                         event.codes.surv = "2",
                         surv.descriptions = "OS",
                         caption = "")
  res2 <- doCoxphGeneric(input.d = lung,
                         var.names = "sex",
                         var.descriptions = "Sex",
                         show.var.detail = TRUE,
                         var.names.surv.time = "time",
                         var.names.surv.status = "status",
                         event.codes.surv = "2",
                         surv.descriptions = "OS",
                         caption = "",
                         var.ref.groups = "2",
                         round.small = TRUE)
  expect_length(res1, 3)
  expect_length(res2, 3)
  expect_identical(res1$result.table[, c(1, 3)], res2$result.table[, c(1, 3)])
  expect_false(res1$result.table[, c(2)] == res2$result.table[, c(2)])
})

test_that("doCoxphGeneric can fit multiple univariable models", {
  res <- doCoxphGeneric(input.d = lung,
                        var.names = c("sex", "ph.ecog"),
                        var.descriptions = c("Sex", "ECOG score"),
                        show.var.detail = TRUE,
                        var.names.surv.time = "time",
                        var.names.surv.status = "status",
                        event.codes.surv = "2",
                        surv.descriptions = "OS",
                        caption = "",
                        var.ref.groups = c("2", "0"))
  expect_equal(nrow(res$result.table), 2)
})

test_that("doCoxphGeneric rounding of small p-values works", {
  res <- doCoxphGeneric(input.d = lung,
                        var.names = c("sex", "ph.ecog"),
                        var.descriptions = c("Sex", "ECOG score"),
                        show.var.detail = TRUE,
                        var.names.surv.time = "time",
                        var.names.surv.status = "status",
                        event.codes.surv = "2",
                        surv.descriptions = "OS",
                        caption = "",
                        round.small = TRUE)
  expect_equal(nrow(res$result.table), 2)
})

test_that("doCoxphGeneric factors use lowest group if no reference specified", {
  lung2 <- lung
  lung2$sex <- as.factor(lung2$sex)
  res1 <- doCoxphGeneric(input.d = lung,
                         var.names = "sex",
                         var.descriptions = "Sex",
                         show.var.detail = TRUE,
                         var.names.surv.time = "time",
                         var.names.surv.status = "status",
                         event.codes.surv = "2",
                         surv.descriptions = "OS",
                         caption = "")
  res2 <- doCoxphGeneric(input.d = lung2,
                         var.names = "sex",
                         var.descriptions = "Sex",
                         show.var.detail = TRUE,
                         var.names.surv.time = "time",
                         var.names.surv.status = "status",
                         event.codes.surv = "2",
                         surv.descriptions = "OS",
                         caption = "")
  expect_identical(res1$result.table, res2$result.table)
})


context("Cox model multivariable")

test_that("doCoxphMultivariable works for univariable case", {
  res <- doCoxphMultivariable(input.d = lung,
                              var.names = "sex",
                              var.descriptions = "Sex",
                              var.names.surv.time = "time",
                              var.names.surv.status = "status",
                              event.codes.surv = "2",
                              surv.descriptions = "OS",
                              caption = "")
  expect_length(res, 4)
})

test_that("doCoxphMultivariable fits multivariable models", {
  res <- doCoxphMultivariable(input.d = lung,
                              var.names = c("sex", "ph.ecog"),
                              var.descriptions = c("Sex", "ECOG score"),
                              show.var.detail = TRUE,
                              var.names.surv.time = "time",
                              var.names.surv.status = "status",
                              event.codes.surv = "2",
                              surv.descriptions = "OS",
                              caption = "",
                              stat.test = "logtest",
                              var.ref.groups = c("2", "0"),
                              round.small = TRUE)
  expect_length(res, 4)
})

test_that("doCoxphMultivariable rounding of small p-values works", {
  res <- doCoxphMultivariable(input.d = lung,
                              var.names = c("sex", "ph.ecog"),
                              var.descriptions = c("Sex", "ECOG score"),
                              show.var.detail = TRUE,
                              var.names.surv.time = "time",
                              var.names.surv.status = "status",
                              event.codes.surv = "2",
                              surv.descriptions = "OS",
                              caption = "",
                              stat.test = "logtest",
                              var.ref.groups = c("2", "0"),
                              round.small = TRUE,
                              round.digits.p.value = 2)
  expect_length(res, 4)
})

test_that(
  "doCoxphMultivariable factors use lowest group if no reference specified", {
    lung2 <- lung
    lung2$sex <- as.factor(lung2$sex)
    res1 <- doCoxphMultivariable(input.d = lung,
                                 var.names = c("sex", "ph.ecog"),
                                 var.descriptions = c("Sex", "ECOG score"),
                                 show.var.detail = TRUE,
                                 var.names.surv.time = "time",
                                 var.names.surv.status = "status",
                                 event.codes.surv = "2",
                                 surv.descriptions = "OS",
                                 caption = "")
    res2 <- doCoxphMultivariable(input.d = lung2,
                                 var.names = c("sex", "ph.ecog"),
                                 var.descriptions = c("Sex", "ECOG score"),
                                 show.var.detail = TRUE,
                                 var.names.surv.time = "time",
                                 var.names.surv.status = "status",
                                 event.codes.surv = "2",
                                 surv.descriptions = "OS",
                                 caption = "")
    expect_identical(res1$result.table, res2$result.table)
  }
)

test_that("doCoxphMultivariable multiple survival outcomes works", {
  res <- doCoxphMultivariable(input.d = lung,
                              var.names = c("sex", "ph.ecog"),
                              var.descriptions = c("Sex", "ECOG score"),
                              show.var.detail = TRUE,
                              var.names.surv.time = c("time", "time2"),
                              var.names.surv.status = c("status", "status2"),
                              event.codes.surv = c("2", "2"),
                              surv.descriptions = c("OS", "OS2"),
                              caption = "")
  expect_equal(nrow(res$result.table), 4)
})


context("Cox model testing interaction")

test_that("doInteractionCox works without an interaction term", {
  res <- doInteractionCox(input.d = lung,
                          var.names = c("sex", "ph.ecog"),
                          var.descriptions = c("Sex", "ECOG score"),
                          var.names.surv.time = c("time", "time2"),
                          var.names.surv.status = c("status", "status2"),
                          event.codes.surv = c("2", "2"),
                          surv.descriptions = c("OS", "OS2"),
                          var.ref.groups = c("2", "0"),
                          caption = "")
  expect_length(res, 4)
})

test_that("doInteractionCox works with an interaction term", {
  lung[["ph.karno:pat.karno"]] <- lung$ph.karno * lung$pat.karno
  lung$ph.karno <- factor(lung$ph.karno)
  lung$pat.karno <- factor(lung$pat.karno)
  res <- suppressWarnings(
    doInteractionCox(input.d = lung,
                     var.names = "ph.karno:pat.karno",
                     var.descriptions = "karnos",
                     var.names.surv.time = c("time", "time2"),
                     var.names.surv.status = c("status", "status2"),
                     event.codes.surv = c("2", "2"),
                     surv.descriptions = c("OS", "OS2"),
                     var.ref.groups = c("2", "0"),
                     caption = "")
  )
  expect_length(res, 4)
})

test_that(
  "doInteractionCox factors use lowest group if no reference specified", {
    lung$sex <- as.factor(lung$sex)
    res <- doInteractionCox(input.d = lung,
                            var.names = c("sex", "ph.ecog"),
                            var.descriptions = c("Sex", "ECOG score"),
                            var.names.surv.time = "time",
                            var.names.surv.status = "status",
                            event.codes.surv = "2",
                            surv.descriptions = "OS",
                            caption = "")
    expect_length(res, 4)
  }
)
