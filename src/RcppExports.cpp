// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// get_cass_vec
std::vector<std::vector<int> > get_cass_vec(std::vector<std::string> c);
RcppExport SEXP _loxcoder_get_cass_vec(SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(get_cass_vec(c));
    return rcpp_result_gen;
END_RCPP
}
// is_valid
std::vector<bool> is_valid(SEXP c);
RcppExport SEXP _loxcoder_is_valid(SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(is_valid(c));
    return rcpp_result_gen;
END_RCPP
}
// decode
Rcpp::S4 decode(std::vector<std::string> r, Rcpp::DataFrame meta, int min_r1_len, int min_r2_len);
RcppExport SEXP _loxcoder_decode(SEXP rSEXP, SEXP metaSEXP, SEXP min_r1_lenSEXP, SEXP min_r2_lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type r(rSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type meta(metaSEXP);
    Rcpp::traits::input_parameter< int >::type min_r1_len(min_r1_lenSEXP);
    Rcpp::traits::input_parameter< int >::type min_r2_len(min_r2_lenSEXP);
    rcpp_result_gen = Rcpp::wrap(decode(r, meta, min_r1_len, min_r2_len));
    return rcpp_result_gen;
END_RCPP
}
// load_origin_files_wrapper
void load_origin_files_wrapper(std::vector<std::string> paths);
RcppExport SEXP _loxcoder_load_origin_files_wrapper(SEXP pathsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type paths(pathsSEXP);
    load_origin_files_wrapper(paths);
    return R_NilValue;
END_RCPP
}
// wrapper_fill_tables
void wrapper_fill_tables();
RcppExport SEXP _loxcoder_wrapper_fill_tables() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    wrapper_fill_tables();
    return R_NilValue;
END_RCPP
}
// pack
std::vector<long long> pack(SEXP c, std::vector<bool> v);
RcppExport SEXP _loxcoder_pack(SEXP cSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type c(cSEXP);
    Rcpp::traits::input_parameter< std::vector<bool> >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(pack(c, v));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_loxcoder_get_cass_vec", (DL_FUNC) &_loxcoder_get_cass_vec, 1},
    {"_loxcoder_is_valid", (DL_FUNC) &_loxcoder_is_valid, 1},
    {"_loxcoder_decode", (DL_FUNC) &_loxcoder_decode, 4},
    {"_loxcoder_load_origin_files_wrapper", (DL_FUNC) &_loxcoder_load_origin_files_wrapper, 1},
    {"_loxcoder_wrapper_fill_tables", (DL_FUNC) &_loxcoder_wrapper_fill_tables, 0},
    {"_loxcoder_pack", (DL_FUNC) &_loxcoder_pack, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_loxcoder(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
