
<!-- README.md is generated from README.Rmd. Please edit that file -->
ageassn
=======

The goal of `ageassn` is to provide a suite of functions to assist in the analysis of the [BrCa Age Associated Biomarker project](https://github.com/BCCRCMO/BrCa_Age_Associated)

Installation
------------

The `ageassn` package is available privately on GitHub. Follow these steps to install directly from GitHub:

1.  Navigate to <https://github.com/settings/tokens> and **Generate new token**
2.  Check off the **repo** scope and then **Generate token**. We have now created a personal access token (PAT) that allows full control of private repositories.
3.  Back at the token home page, write down the PAT you just generated and save the following in a file named `.Renviron`, making sure the last line is a newline:

    ``` r
    GITHUB_PAT="YOUR_PERSONAL_ACCESS_TOKEN"
    ```

4.  Put the `.Renviron` file in the project root directory where you want to install the package (e.g. project directory for the age paper or even this package!) and add it to your `.gitignore` so you protect your PAT
5.  Now we can install the package just like from a public repository:

``` r
# install.packages("devtools")
devtools::install_github("BCCRCMO/ageassn")
```
