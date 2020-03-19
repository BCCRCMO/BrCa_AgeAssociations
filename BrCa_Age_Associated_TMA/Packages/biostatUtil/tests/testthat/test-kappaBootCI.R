
context("Kappa bootstrapped confidence interval")

s <- 12
set.seed(s)
a <- rbinom(n = 100, size = 1, prob = 0.3)
b <- rbinom(n = 100, size = 1, prob = 0.7)
nb <- 100

test_that("different methods of calculating kappa", {
  expect_error(kappaBootCI(a, b, s, nb, method = "cohen"), NA)
  expect_error(kappaBootCI(a, b, s, nb, method = "weighted"), NA)
  expect_error(kappaBootCI(a, b, s, nb, method = "fleiss"), NA)
  expect_error(kappaBootCI(a, b, s, nb, method = "krippendorff"), NA)
})

test_that("confidence interval labels are correct", {
  cl <- 0.90
  kb <- kappaBootCI(a, b, s, nb, cl, method = "cohen")
  expect_equal(names(kb)[-1],
               c(paste0((1 - cl) / 2 * 100, "%"),
                 paste0((1 - (1 - cl) / 2) * 100, "%")))
})
