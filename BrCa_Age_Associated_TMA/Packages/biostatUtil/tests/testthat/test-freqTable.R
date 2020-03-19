
context("Generating frequency tables and histograms")

# Create vector of randomly reordered alphabet with various frequencies
# for each letter
set.seed(123)
n <- sample(10, length(letters), replace = TRUE)
x <- sample(rep.int(letters, times = n))

test_that("Score label and description label match up", {
  ftdf <- freqTable(x, levels = letters, missing = "x", description = letters,
                    round = 3, plot = TRUE)
  idxp <- ftdf[, "Score"] != "Total" & ftdf[, "Class"] != "Total"
  expect_length(ftdf, 7)
  expect_equal(nrow(ftdf), length(letters) + 2)
  expect_identical(ftdf[idxp, "Score"], ftdf[idxp, "Description"])
  if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
})

test_that("Different layout if no missing cases", {
  ftdf2 <- freqTable(x, levels = letters, description = letters,
                     round = 3, plot = TRUE)
  expect_length(ftdf2, 5)
  expect_equal(nrow(ftdf2), length(letters) + 1)
})
