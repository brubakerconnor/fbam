// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// l2assign
IntegerVector l2assign(NumericMatrix spec, NumericMatrix p);
RcppExport SEXP _fbam_l2assign(SEXP specSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type spec(specSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(l2assign(spec, p));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fbam_l2assign", (DL_FUNC) &_fbam_l2assign, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_fbam(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}