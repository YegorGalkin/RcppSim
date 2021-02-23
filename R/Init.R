.onLoad <- function(libname, pkgname) {
  require(Rcpp)
  require(devtools)
  loadModule("poisson_1d_module",TRUE)
  loadModule("poisson_2d_module",TRUE)
  loadModule("poisson_3d_module",TRUE)
  loadModule("poisson_n_species_module", TRUE)
}
