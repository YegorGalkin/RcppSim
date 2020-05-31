.onLoad <- function(libname, pkgname) {
  require(Rcpp)
  require(devtools)
  loadModule("poisson_1d_module",TRUE)
  loadModule("poisson_2d_module",TRUE)
}
