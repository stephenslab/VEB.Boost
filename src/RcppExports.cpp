// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "VEB.Boost_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// make_stumps_matrix_cpp
XPtr<stumpsmatrix::StumpsMatrix> make_stumps_matrix_cpp(arma::mat X, arma::uvec include_linear, arma::uvec include_stumps, std::vector<arma::vec> cuts, unsigned int ncores);
RcppExport SEXP _VEB_Boost_make_stumps_matrix_cpp(SEXP XSEXP, SEXP include_linearSEXP, SEXP include_stumpsSEXP, SEXP cutsSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type include_linear(include_linearSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type include_stumps(include_stumpsSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::vec> >::type cuts(cutsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(make_stumps_matrix_cpp(X, include_linear, include_stumps, cuts, ncores));
    return rcpp_result_gen;
END_RCPP
}
// make_stumps_matrix_sp_cpp
XPtr<stumpsmatrix::StumpsMatrix> make_stumps_matrix_sp_cpp(arma::sp_mat X, arma::uvec include_linear, arma::uvec include_stumps, std::vector<arma::vec> cuts, unsigned int ncores);
RcppExport SEXP _VEB_Boost_make_stumps_matrix_sp_cpp(SEXP XSEXP, SEXP include_linearSEXP, SEXP include_stumpsSEXP, SEXP cutsSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type include_linear(include_linearSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type include_stumps(include_stumpsSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::vec> >::type cuts(cutsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(make_stumps_matrix_sp_cpp(X, include_linear, include_stumps, cuts, ncores));
    return rcpp_result_gen;
END_RCPP
}
// weighted_SER_cpp
List weighted_SER_cpp(XPtr<stumpsmatrix::StumpsMatrix> xp, arma::vec& Y, arma::vec& sigma2, Nullable<List> init, double max_lV, double lin_prior_prob);
RcppExport SEXP _VEB_Boost_weighted_SER_cpp(SEXP xpSEXP, SEXP YSEXP, SEXP sigma2SEXP, SEXP initSEXP, SEXP max_lVSEXP, SEXP lin_prior_probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<stumpsmatrix::StumpsMatrix> >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type init(initSEXP);
    Rcpp::traits::input_parameter< double >::type max_lV(max_lVSEXP);
    Rcpp::traits::input_parameter< double >::type lin_prior_prob(lin_prior_probSEXP);
    rcpp_result_gen = Rcpp::wrap(weighted_SER_cpp(xp, Y, sigma2, init, max_lV, lin_prior_prob));
    return rcpp_result_gen;
END_RCPP
}
// predFnSusieStumps_cpp
arma::vec predFnSusieStumps_cpp(XPtr<stumpsmatrix::StumpsMatrix> xp_new, List currentFit, int moment);
RcppExport SEXP _VEB_Boost_predFnSusieStumps_cpp(SEXP xp_newSEXP, SEXP currentFitSEXP, SEXP momentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<stumpsmatrix::StumpsMatrix> >::type xp_new(xp_newSEXP);
    Rcpp::traits::input_parameter< List >::type currentFit(currentFitSEXP);
    Rcpp::traits::input_parameter< int >::type moment(momentSEXP);
    rcpp_result_gen = Rcpp::wrap(predFnSusieStumps_cpp(xp_new, currentFit, moment));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_VEB_Boost_make_stumps_matrix_cpp", (DL_FUNC) &_VEB_Boost_make_stumps_matrix_cpp, 5},
    {"_VEB_Boost_make_stumps_matrix_sp_cpp", (DL_FUNC) &_VEB_Boost_make_stumps_matrix_sp_cpp, 5},
    {"_VEB_Boost_weighted_SER_cpp", (DL_FUNC) &_VEB_Boost_weighted_SER_cpp, 6},
    {"_VEB_Boost_predFnSusieStumps_cpp", (DL_FUNC) &_VEB_Boost_predFnSusieStumps_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_VEB_Boost(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
