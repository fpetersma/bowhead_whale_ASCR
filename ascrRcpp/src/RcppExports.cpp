// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// uncondNllRcpp
double uncondNllRcpp(List dat, List par);
RcppExport SEXP _ascrRcpp_uncondNllRcpp(SEXP datSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type dat(datSEXP);
    Rcpp::traits::input_parameter< List >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(uncondNllRcpp(dat, par));
    return rcpp_result_gen;
END_RCPP
}
// singleNllRcpp
double singleNllRcpp(List dat, List par);
RcppExport SEXP _ascrRcpp_singleNllRcpp(SEXP datSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type dat(datSEXP);
    Rcpp::traits::input_parameter< List >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(singleNllRcpp(dat, par));
    return rcpp_result_gen;
END_RCPP
}
// uncondNllRcppFixedSL
double uncondNllRcppFixedSL(List dat, List par);
RcppExport SEXP _ascrRcpp_uncondNllRcppFixedSL(SEXP datSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type dat(datSEXP);
    Rcpp::traits::input_parameter< List >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(uncondNllRcppFixedSL(dat, par));
    return rcpp_result_gen;
END_RCPP
}
// singleNllRcppFixedSL
double singleNllRcppFixedSL(List dat, List par);
RcppExport SEXP _ascrRcpp_singleNllRcppFixedSL(SEXP datSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type dat(datSEXP);
    Rcpp::traits::input_parameter< List >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(singleNllRcppFixedSL(dat, par));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ascrRcpp_uncondNllRcpp", (DL_FUNC) &_ascrRcpp_uncondNllRcpp, 2},
    {"_ascrRcpp_singleNllRcpp", (DL_FUNC) &_ascrRcpp_singleNllRcpp, 2},
    {"_ascrRcpp_uncondNllRcppFixedSL", (DL_FUNC) &_ascrRcpp_uncondNllRcppFixedSL, 2},
    {"_ascrRcpp_singleNllRcppFixedSL", (DL_FUNC) &_ascrRcpp_singleNllRcppFixedSL, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ascrRcpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
