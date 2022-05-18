// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// d_2DSD
NumericVector d_2DSD(NumericVector rts, NumericVector params, double precision, int boundary, bool stop_on_error);
RcppExport SEXP _dynConfiR_d_2DSD(SEXP rtsSEXP, SEXP paramsSEXP, SEXP precisionSEXP, SEXP boundarySEXP, SEXP stop_on_errorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rts(rtsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< double >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< int >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< bool >::type stop_on_error(stop_on_errorSEXP);
    rcpp_result_gen = Rcpp::wrap(d_2DSD(rts, params, precision, boundary, stop_on_error));
    return rcpp_result_gen;
END_RCPP
}
// d_WEVmu
NumericVector d_WEVmu(NumericVector rts, NumericVector params, double precision, int boundary, bool stop_on_error);
RcppExport SEXP _dynConfiR_d_WEVmu(SEXP rtsSEXP, SEXP paramsSEXP, SEXP precisionSEXP, SEXP boundarySEXP, SEXP stop_on_errorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rts(rtsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< double >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< int >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< bool >::type stop_on_error(stop_on_errorSEXP);
    rcpp_result_gen = Rcpp::wrap(d_WEVmu(rts, params, precision, boundary, stop_on_error));
    return rcpp_result_gen;
END_RCPP
}
// d_IRM
NumericVector d_IRM(NumericVector rts, NumericVector params, int win, double step_width);
RcppExport SEXP _dynConfiR_d_IRM(SEXP rtsSEXP, SEXP paramsSEXP, SEXP winSEXP, SEXP step_widthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rts(rtsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< double >::type step_width(step_widthSEXP);
    rcpp_result_gen = Rcpp::wrap(d_IRM(rts, params, win, step_width));
    return rcpp_result_gen;
END_RCPP
}
// d_PCRM
NumericVector d_PCRM(NumericVector rts, NumericVector params, int win, double step_width);
RcppExport SEXP _dynConfiR_d_PCRM(SEXP rtsSEXP, SEXP paramsSEXP, SEXP winSEXP, SEXP step_widthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rts(rtsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< double >::type step_width(step_widthSEXP);
    rcpp_result_gen = Rcpp::wrap(d_PCRM(rts, params, win, step_width));
    return rcpp_result_gen;
END_RCPP
}
// dd_IRM
NumericVector dd_IRM(NumericVector rts, NumericVector xj, NumericVector params, int win, int method);
RcppExport SEXP _dynConfiR_dd_IRM(SEXP rtsSEXP, SEXP xjSEXP, SEXP paramsSEXP, SEXP winSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rts(rtsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xj(xjSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(dd_IRM(rts, xj, params, win, method));
    return rcpp_result_gen;
END_RCPP
}
// dd_PCRM
NumericVector dd_PCRM(NumericVector rts, NumericVector xj, NumericVector params, int win);
RcppExport SEXP _dynConfiR_dd_PCRM(SEXP rtsSEXP, SEXP xjSEXP, SEXP paramsSEXP, SEXP winSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rts(rtsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xj(xjSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    rcpp_result_gen = Rcpp::wrap(dd_PCRM(rts, xj, params, win));
    return rcpp_result_gen;
END_RCPP
}
// r_RM
NumericVector r_RM(int n, NumericVector params, double rho, double delta, double maxT);
RcppExport SEXP _dynConfiR_r_RM(SEXP nSEXP, SEXP paramsSEXP, SEXP rhoSEXP, SEXP deltaSEXP, SEXP maxTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type maxT(maxTSEXP);
    rcpp_result_gen = Rcpp::wrap(r_RM(n, params, rho, delta, maxT));
    return rcpp_result_gen;
END_RCPP
}
// r_WEV
NumericVector r_WEV(int n, NumericVector params, int model, double delta, double maxT, bool stop_on_error);
RcppExport SEXP _dynConfiR_r_WEV(SEXP nSEXP, SEXP paramsSEXP, SEXP modelSEXP, SEXP deltaSEXP, SEXP maxTSEXP, SEXP stop_on_errorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type maxT(maxTSEXP);
    Rcpp::traits::input_parameter< bool >::type stop_on_error(stop_on_errorSEXP);
    rcpp_result_gen = Rcpp::wrap(r_WEV(n, params, model, delta, maxT, stop_on_error));
    return rcpp_result_gen;
END_RCPP
}
// r_RM_Kiani
NumericVector r_RM_Kiani(int n, NumericVector params, double rho, double Bl, double delta, double maxT);
RcppExport SEXP _dynConfiR_r_RM_Kiani(SEXP nSEXP, SEXP paramsSEXP, SEXP rhoSEXP, SEXP BlSEXP, SEXP deltaSEXP, SEXP maxTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type Bl(BlSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type maxT(maxTSEXP);
    rcpp_result_gen = Rcpp::wrap(r_RM_Kiani(n, params, rho, Bl, delta, maxT));
    return rcpp_result_gen;
END_RCPP
}
// r_LCA
NumericVector r_LCA(int n, NumericVector params, double delta, double maxT);
RcppExport SEXP _dynConfiR_r_LCA(SEXP nSEXP, SEXP paramsSEXP, SEXP deltaSEXP, SEXP maxTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type maxT(maxTSEXP);
    rcpp_result_gen = Rcpp::wrap(r_LCA(n, params, delta, maxT));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dynConfiR_d_2DSD", (DL_FUNC) &_dynConfiR_d_2DSD, 5},
    {"_dynConfiR_d_WEVmu", (DL_FUNC) &_dynConfiR_d_WEVmu, 5},
    {"_dynConfiR_d_IRM", (DL_FUNC) &_dynConfiR_d_IRM, 4},
    {"_dynConfiR_d_PCRM", (DL_FUNC) &_dynConfiR_d_PCRM, 4},
    {"_dynConfiR_dd_IRM", (DL_FUNC) &_dynConfiR_dd_IRM, 5},
    {"_dynConfiR_dd_PCRM", (DL_FUNC) &_dynConfiR_dd_PCRM, 4},
    {"_dynConfiR_r_RM", (DL_FUNC) &_dynConfiR_r_RM, 5},
    {"_dynConfiR_r_WEV", (DL_FUNC) &_dynConfiR_r_WEV, 6},
    {"_dynConfiR_r_RM_Kiani", (DL_FUNC) &_dynConfiR_r_RM_Kiani, 6},
    {"_dynConfiR_r_LCA", (DL_FUNC) &_dynConfiR_r_LCA, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_dynConfiR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
