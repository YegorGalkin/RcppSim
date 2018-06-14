# RcppSim
Rcpp interface for c++ based poisson simulator

Requires Rtools installed

R packages: Rcpp, BH, devtools

Install with:

devtools::install_github("YegorGalkin/RcppSim",quick=TRUE,local=FALSE)


How to get interactive version running:

Install docker

Run in console, WORKING_DIR is a diretory to be mounted in jupyter home dir:

docker run -p 8888:8888 -v WORKING_DIR:/home/jovyan/work jupyter/datascience-notebook

Copy one time token

Login on http://localhost:8888

Paste one time token

Run in R kernel:

install.packages("Rcpp")

require(devtools)

require(Rcpp)

require(BH)

require(testthat)

options(unzip = "internal")

devtools::install_github("YegorGalkin/RcppSim",quick=TRUE,local=FALSE)

Thats it, simulator is good to go!
