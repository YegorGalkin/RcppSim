#' @useDynLib MathBioSim, .registration = TRUE
#' @export poisson_1d
#' @export poisson_2d
#' @export poisson_3d
#' @import Rcpp
NULL

Rcpp::loadModule("poisson_1d_module", TRUE)
Rcpp::loadModule("poisson_2d_module", TRUE)
Rcpp::loadModule("poisson_3d_module", TRUE)
