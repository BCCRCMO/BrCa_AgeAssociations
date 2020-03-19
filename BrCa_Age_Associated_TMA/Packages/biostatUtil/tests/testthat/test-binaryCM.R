
context("Binary class confusion matrices")

set.seed(547)
n <- 80
tol <- .0002
x <- rbinom(n, size = 1, prob = 0.6) # Ref
xf <- factor(x, labels = c("positive", "negative"))

y <- rbinom(n, size = 1, prob = 0.4) # Pred
yf <- factor(y, labels = c("negative", "positive"))

expect_binaryCM_success <- function(x, y, pcond) {
  bcm <- binaryCM(x, y, pcond = pcond)
  z <- apply(rmultinom(n, size = 1, prob = c(0.2, 0.4, 0.4)), 2,
             function(x) which(x == 1))

  test_that("categorical variable > 2 levels will result in an error", {
    expect_error(binaryCM(x, z))
  })

  test_that("positive condition is in (1,1) position of confusion matrix", {
    expect_equal(sum(x == pcond & y == pcond), bcm$CM[1, 1])
  })

  test_that("the output of the table corresponds to the individual elements", {
    expect_true(all(bcm$table["Accuracy", ] == bcm$Accuracy))
    expect_true(all(bcm$table["Sensitivity", ] == bcm$Sensitivity))
    expect_true(all(bcm$table["Specificity", ] == bcm$Specificity))
    expect_true(all(bcm$table["PPV", ] == bcm$PPV))
    expect_true(all(bcm$table["NPV", ] == bcm$NPV))
    expect_true(all(bcm$table["kappa", ] == bcm$kappa))
  })

  test_that("the point estimate is included in the 95 bootstrap interval", {
    expect_true(all(bcm$table[, "Point Estimate"] > bcm$table[, "Lower CI"]) &
                  all(bcm$table[, "Point Estimate"] < bcm$table[, "Upper CI"]))
  })

  test_that("accuracy is sum of diagonal over total", {
    expect_equal(
      bcm$Accuracy[, "PointEst"],
      (sum(x == pcond & y == pcond) + sum(x != pcond & y != pcond)) / length(x)
    )
  })

  test_that("sens is sum of true positives / sum of true condition positive", {
    expect_equal(
      bcm$Sensitivity[, "PointEst"],
      sum(x == pcond & y == pcond) / sum(x == pcond), tolerance = tol
    )
  })

  test_that("spec is sum of true negatives / sum of true condition negative", {
    expect_equal(
      bcm$Specificity[, "PointEst"],
      sum(x != pcond & y != pcond) / sum(x != pcond), tolerance = tol
    )
  })

  test_that("PPV is sum of true positives / predicted positive", {
    expect_equal(
      bcm$PPV[, "PointEst"],
      sum(x == pcond & y == pcond) / sum(y == pcond), tolerance = tol
    )
  })

  test_that("NPV is sum of true negatives / predicted negative", {
    expect_equal(
      bcm$NPV[, "PointEst"],
      sum(x != pcond & y != pcond) / sum(y != pcond), tolerance = tol
    )
  })
}

test_that("binaryCM works for numeric values (0/1)", {
  expect_error(expect_binaryCM_success(x, y, pcond = 1), NA)
})

test_that("binaryCM works with factors and labels", {
  expect_error(expect_binaryCM_success(xf, yf, pcond = "positive"), NA)
})

test_that("outputs can print to console", {
  expect_output(binaryCM(x, y, verbose = TRUE))
})
