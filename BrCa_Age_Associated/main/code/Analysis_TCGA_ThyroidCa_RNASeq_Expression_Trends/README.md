# TCGA Thyroid RNA-Seq Expression Trends

Analyses are separated into [AllCases](https://github.com/BCCRCMO/BrCa_Age_Associated/tree/master/main/code/Analysis_TCGA_ThyroidCa_RNASeq_Expression_Trends/AllCases), [Females](https://github.com/BCCRCMO/BrCa_Age_Associated/tree/master/main/code/Analysis_TCGA_ThyroidCa_RNASeq_Expression_Trends/Females), and [Males](https://github.com/BCCRCMO/BrCa_Age_Associated/tree/master/main/code/Analysis_TCGA_ThyroidCa_RNASeq_Expression_Trends/Males). In each of these directories, there are subdirectories `analysis`, `figures`, and `tables.`

There are four biological/statistical significance levels considered:

- FDR = 0.05, FC = 1.25
- FDR = 0.01, FC = 1.25
- FDR = 0.01, FC = 2
- FDR = 0.02, FC = 4

The contents are:

- `analysis`: Most analytical results are produced by `TCGA_{organ}_Study*.R`. Patient characteristics are summarized in `*PtChars.R`. There is also a Fisher Exact Test report for both ER binding definitions.
- `figures`: For each level, probe level, volcano, and manhattan plots. In `single_figures`, there are individual probe level plots. Version 1 volcano plots (`AgeRelated_Volcano*.pdf`) labels on the original ER binding status, and version 2 volcano plots (`AgeRelated_Volcano*v02.pdf`) labels on the NKI ER binding Tier 1/2 definitions.
- `tables`: For each level, `aroutdf*.csv` contains the regression statistics on all probes, `AgeDependent_aroutdf*.csv` filters `aroutdf*.csv` for age-dependent probes, and `BHadj_and_AgeDependent_aroutdf*.csv` filters for BH-adjusted p-values *and* age-dependence.