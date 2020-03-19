context("Independence Tests")

library(descr)

test_that("Output works with defaults", {
  set.seed(1108)
  A <- rbinom(100, 3, 0.2)
  B <- rbinom(100, 4, 0.8)
  ct <- CrossTable(A, B, prop.chisq = FALSE, chisq = FALSE, fisher = TRUE)
  expect_error(indepTests(ct), NA)
})

test_that("Continuity Correction only when Pearson is not NA", {
  set.seed(1108)
  A <- rbinom(100, 1, 0.2)
  B <- rbinom(100, 1, 0.8)
  expect_warning(indepTests(CrossTable(A, B)))

  ct2 <- CrossTable(matrix(c(20, 25, 15, 30), nrow = 2))
  df <- indepTests(ct2)
  expect_false(anyNA(df[1:2, -1]))
})

test_that("NA are returned when statistics can't be computed", {
  ct <- CrossTable(matrix(c(0, 0, 0, 0), nrow = 2), prop.chisq = FALSE)
  df <- indepTests(ct)
  expect_true(all(is.na(df[1:5, -1])))
})
