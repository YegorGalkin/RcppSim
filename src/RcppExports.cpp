// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;


RcppExport SEXP _rcpp_module_boot_poisson_1d_module();
RcppExport SEXP _rcpp_module_boot_poisson_2d_module();
RcppExport SEXP _rcpp_module_boot_poisson_3d_module();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_poisson_1d_module", (DL_FUNC) &_rcpp_module_boot_poisson_1d_module, 0},
    {"_rcpp_module_boot_poisson_2d_module", (DL_FUNC) &_rcpp_module_boot_poisson_2d_module, 0},
    {"_rcpp_module_boot_poisson_3d_module", (DL_FUNC) &_rcpp_module_boot_poisson_3d_module, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_MathBioSim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}