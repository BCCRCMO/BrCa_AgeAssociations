# biostatUtil installation

When updating R to a new version, there are several administrative tasks to perform to ensure the package can still be installed.

1. Run `devtools::install()` to install all dependencies. These are all packages listed in the `DESCRIPTION` `Imports` field
2. Install any remaining dependencies from Bioconductor
3. Run `devtools::check()` to see which packages to install from `DESCRIPTION` `Suggests` field
4. Open Terminal and run `R CMD javareconf` to reconfigure Java paths and other configurations. Note that Java v1.8.0 is required for java dependencies to successfully load when `biostatUtil` is attached.
