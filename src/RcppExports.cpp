// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// decode
Rcpp::S4 decode(std::vector<std::string> r, Rcpp::DataFrame meta);
RcppExport SEXP _loxcoder_decode(SEXP rSEXP, SEXP metaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type r(rSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type meta(metaSEXP);
    rcpp_result_gen = Rcpp::wrap(decode(r, meta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_loxcoder_decode", (DL_FUNC) &_loxcoder_decode, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_loxcoder(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
