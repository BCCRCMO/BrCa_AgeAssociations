
context("Print streamlined logistf output")

library(logistf)
data(sex2)
fit <- logistf(case ~ age + oc + vic + vicl + vis + dia, data = sex2,
               alpha = 0.1)

test_that("streamlined summary is a matrix", {
  expect_is(logistfPrint(fit), "matrix")
})
