.onLoad <- function(libname, pkgname) {
  require(Rcpp)
  require(devtools)
  loadModule("poisson_1d_module",TRUE)
  loadModule("poisson_1d_n_species_module",TRUE)
}
