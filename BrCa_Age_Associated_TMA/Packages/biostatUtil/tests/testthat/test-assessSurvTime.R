
context("Assess survival time")

set.seed(3)
starts <- seq(as.Date("2000/1/1"), as.Date("2003/1/1"), by = "quarter")
ends <-  seq(as.Date("2003/1/1"), as.Date("2006/1/1"), by = "quarter")
statuses <- sample(0:1, 13, replace = TRUE)
x1 <- assessSurvTime(starts, ends, statuses)
statuses[9:12] <- NA
x2 <- assessSurvTime(starts, ends, statuses)

test_that("output is a list of 5 real-valued elements", {
  expect_is(x1, "list")
  expect_length(x1, 5)
})

test_that("missing cases are handled", {
  expect_false(isTRUE(all.equal(x1, x2)))
})
