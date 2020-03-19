
context("Mass Spectrometry functions")

data(psm, pro)
dat1 <- ms_process(psm = psm, protein = pro, treatment = c("N-FLAG", "UNTAGGED"),
                   control = "EMPTY", condition = "MRC",
                   sample.id = c("CLK1_N-FLAG_r1", "CLK1_UNTAGGED_r1", "FLAG_EMPTY_r1",
                                 "CLK1_N-FLAG_r2", "CLK1_UNTAGGED_r2", "FLAG_EMPTY_r2"))
dat2 <- ms_process(psm = psm, protein = pro, treatment = c("EMPTY", "N-FLAG", "UNTAGGED"),
                   condition = "MRC",
                   sample.id = c("CLK1_N-FLAG_r1", "CLK1_UNTAGGED_r1", "FLAG_EMPTY_r1",
                                 "CLK1_N-FLAG_r2", "CLK1_UNTAGGED_r2", "FLAG_EMPTY_r2"))
dat3 <- ms_process(psm = psm, protein = pro, treatment = c("1", "2", "3"),
                   condition = "VC")

test_that("ms_process works with/without controls and sample ids", {
  expect_error(dat1, NA)
  expect_error(dat2, NA)
  expect_error(dat3, NA)
})

test_that("ms_boxplot shows boxplots of vsn values across samples", {
  bp <- ms_boxplot(dat1)
  expect_length(bp, 3)
  expect_equal(unique(unlist(lapply(bp, class))), c("gg", "ggplot"))
})

# Vector of groups
g <- c("EMPTY", "N-FLAG", "UNTAGGED")

test_that("ms_mean_var shows mean variance relationship", {
  mv <- ms_mean_var(dat1, g = g,
                    title = c("No Flag", "CLK1 N-Flagged", "CLK1 Untagged"))
  expect_length(mv, length(g) + 1)
  expect_equal(unique(unlist(lapply(mv, class))), c("gg", "ggplot"))
})

test_that("ms_mean_var takes title from g if none provided", {
  mv_no_title <- ms_mean_var(dat1, g = g)
  expect_equal(purrr::map_chr(mv_no_title, c("labels", "title")),
               c("vsn(Raw data values)", paste("vsn", g)))
})

all_vars <-
  c("Gene",
    "Omnibus Treatment F statistic",
    "Omnibus Treatment F statistic numerator degrees of freedom",
    "Omnibus Treatment F statistic denominator degrees of freedom",
    "Omnibus Treatment p-value",
    "Benjamini Hochberg adjusted omnibus p-value",
    "Fitted Mean of Empty group",
    "Fitted Mean of N-Flagged group",
    "Fitted Mean of Untagged group",
    "Standard Error of full model",

    "N-Flagged Average treatment effect t statistic value",
    "N-Flagged Average treatment effect t statistic Wald-based p-value",
    "N-Flagged Benjamini Hochberg adjusted Wald p-value for average treatment effect",
    "N-Flagged log2(average treatment effect difference from control) lower 95% confidence limit",
    "N-Flagged log2(average treatment effect difference from control)",
    "N-Flagged log2(average treatment effect difference from control) upper 95% confidence limit",
    "N-Flagged Fold change lower 95% confidence limit",
    "N-Flagged Fold change for average treatment effect",
    "N-Flagged Fold change upper 95% confidence limit",
    "N-Flagged Absolute fold change direction",
    "N-Flagged Absolute fold change lower 95% confidence interval limit",
    "N-Flagged Absolute fold change (average treatment effect across peptides)",
    "N-Flagged Absolute fold change upper 95% confidence interval limit",

    "Untagged Average treatment effect t statistic value",
    "Untagged Average treatment effect t statistic Wald-based p-value",
    "Untagged Benjamini Hochberg adjusted Wald p-value for average treatment effect",
    "Untagged log2(average treatment effect difference from control) lower 95% confidence limit",
    "Untagged log2(average treatment effect difference from control)",
    "Untagged log2(average treatment effect difference from control) upper 95% confidence limit",
    "Untagged Fold change lower 95% confidence limit",
    "Untagged Fold change for average treatment effect",
    "Untagged Fold change upper 95% confidence limit",
    "Untagged Absolute fold change direction",
    "Untagged Absolute fold change lower 95% confidence interval limit",
    "Untagged Absolute fold change (average treatment effect across peptides)",
    "Untagged Absolute fold change upper 95% confidence interval limit",

    "Accession",
    "Sequence",
    "Annotated.Sequence",
    "Descriptions",
    "Modifications",
    "Reporter.Quan.Result.ID")
info_vars <- c("Gene", "Accession", "Sequence", "Annotated.Sequence",
               "Descriptions", "Modifications", "Reporter.Quan.Result.ID")
resMat <- suppressWarnings(
  ms_summarize(dat1, g = c("FLAG_EMPTY", "CLK1_N-FLAG", "CLK1_UNTAGGED"),
               level = "Gene", col.names = all_vars, info.vars = info_vars))

test_that("ms_summarize works", {
  expect_error(resMat, NA)
  expect_length(resMat, 14 * length(g))
})

test_that("ms_top works", {
  topMat <- ms_top(resMat, level = "Gene")
  expect_error(topMat, NA)
})
