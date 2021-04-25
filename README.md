# Poisson point process simulation of spatial population dynamics

This repository contains source code (R/Rcpp) required to run spatial population dynamics simualation.

# Installation
### Windows
  - Install R ([R Installer](https://cran.r-project.org/bin/windows/base/)) - latest version used is 4.0.5
  - Install Rtools ([Rtools Installer](https://cran.r-project.org/bin/windows/Rtools/))
  - Install required R packages from CRAN (devtools, Rcpp, BH)
  - Install this source package using devtools
```R
devtools::install_github("YegorGalkin/RcppSim")
```
# Usage
See examples for simulator usage in "examples" folder
# Changelog
## April 2021
 - Removed population cap parameter on simulations
 - Added realtime limit parameter on simulations
 - Birth kernel is now passed as quantile function for radius distribution
 - Added R function wrappers for initializing and running simulators
 - Added examples for running 1d simulations and parallel execution
 - Added some parameter value checks and default values for model
 - Added experimental function for batch running simulations
 - Standartized interfaces between 1d,2d and 3d simulations
 - Added documentation for functions via roxygen2
 - Fixed package initialization and function registration (Rcpp)
