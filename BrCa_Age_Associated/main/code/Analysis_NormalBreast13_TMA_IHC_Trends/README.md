# NormalBreast13 TMA IHC Trends

The Normal Breast 13 TMA IHC Trends analyses were carried out on the remote `ssh://git@10.9.18.229/Users/git/TMA`. Summary [Plots](https://github.com/BCCRCMO/BrCa_Age_Associated/tree/master/main/code/Analysis_NormalBreast13_TMA_IHC_Trends/Plots) and [Tables](https://github.com/BCCRCMO/BrCa_Age_Associated/tree/master/main/code/Analysis_NormalBreast13_TMA_IHC_Trends/Tables) were copied over to this repository.

IHC trend analysis used the Friedman supersmoother and its confidence interval was generated using 4 different significance/simulation levels:

- 95% CI at 2000 simulations
- 99% CI at 10000 simulations
- 99.9% CI at 40000 simulations
- 99.99% CI at 400000 simulations

## Plots

The analysis was performed on WholeSeries of the full NormalBreast13 TMA data. Figures with results for all 12 biomarkers are in the top level of [Plots](https://github.com/BCCRCMO/BrCa_Age_Associated/tree/master/main/code/Analysis_NormalBreast13_TMA_IHC_Trends/Plots), but individual biomarker figures can be found in the [subdirectory](https://github.com/BCCRCMO/BrCa_Age_Associated/tree/master/main/code/Analysis_NormalBreast13_TMA_IHC_Trends/Plots/NormalBreast13_IHC_Trends_Individual).

## Tables

There is a "logical" [significance table](https://github.com/BCCRCMO/BrCa_Age_Associated/blob/master/main/code/Analysis_NormalBreast13_TMA_IHC_Trends/Tables/NormalBreast13_IHC_Significance.csv). Each row represents a gene, and a column for each of the 4 significance levels. A `TRUE` indicates that the gene at that particular level was age-dependent, and `FALSE` otherwise.

A [correlation table](https://github.com/BCCRCMO/BrCa_Age_Associated/blob/master/main/code/Analysis_NormalBreast13_TMA_IHC_Trends/Tables/NormalBreast13_TMA_Correlations.csv) summarizes the proportion of the assay range covered on average by the fitted supersmoother values. Specifically, they are calculated as

`(smoother value at max age - smoother value at min age) / (assay max value - assay min value)`
