#include "MatrixClass.h"
#include <RcppArmadillo.h>
#include <memory>

using namespace arma;
using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
XPtr<stumpsmatrix::StumpsMatrix> make_stumps_matrix_cpp(arma::mat& X, Rcpp::Nullable<Rcpp::IntegerVector> include_linear, Rcpp::Nullable<Rcpp::IntegerVector> include_stumps, unsigned int num_cuts, bool use_quants, unsigned int scale_X, unsigned int ncores) {
  stumpsmatrix::StumpsMatrix *X_stumps = new stumpsmatrix::StumpsMatrix(X, include_linear, include_stumps, num_cuts, use_quants, scale_X, ncores);
  XPtr<stumpsmatrix::StumpsMatrix> X_stumps_p(X_stumps, true);
  return(X_stumps_p);
}

// [[Rcpp::export]]
XPtr<stumpsmatrix::StumpsMatrix> make_stumps_matrix_sp_cpp(arma::sp_mat& X, Rcpp::Nullable<Rcpp::IntegerVector> include_linear, Rcpp::Nullable<Rcpp::IntegerVector> include_stumps, unsigned int num_cuts, bool use_quants, unsigned int scale_X, unsigned int ncores) {
  stumpsmatrix::StumpsMatrix *X_stumps = new stumpsmatrix::StumpsMatrix(X, include_linear, include_stumps, num_cuts, use_quants, scale_X, ncores);
  XPtr<stumpsmatrix::StumpsMatrix> X_stumps_p(X_stumps, true);
  return(X_stumps_p);
}

// [[Rcpp::export]]
XPtr<stumpsmatrix::StumpsMatrix> make_stumps_test_matrix_cpp(arma::mat& X, XPtr<stumpsmatrix::StumpsMatrix> SM_in) {
    stumpsmatrix::StumpsMatrix* X_stumps = new stumpsmatrix::StumpsMatrix(X, SM_in);
    XPtr<stumpsmatrix::StumpsMatrix> X_stumps_p(X_stumps, true);
    return(X_stumps_p);
}

// [[Rcpp::export]]
XPtr<stumpsmatrix::StumpsMatrix> make_stumps_test_matrix_sp_cpp(arma::sp_mat& X, XPtr<stumpsmatrix::StumpsMatrix> SM_in) {
    stumpsmatrix::StumpsMatrix* X_stumps = new stumpsmatrix::StumpsMatrix(X, SM_in);
    XPtr<stumpsmatrix::StumpsMatrix> X_stumps_p(X_stumps, true);
    return(X_stumps_p);
}

// compute KL-divergence
double calc_KL_cpp(vec& mu, vec& alpha, vec& sigma2, vec& prior_weights, double V) {
  vec KL_div = alpha % (log((alpha / prior_weights) % sqrt(V / sigma2)) - .5 + ((sigma2 + pow(mu, 2)) / (2 * V)));
  KL_div.elem( find(alpha == 0.0) ).zeros();
  // KL_div.elem( find(approx_equal(alpha, 0.0, "absdiff", 1e-12)) ) = 0.0;
  // KL_div.elem( find(find_nonfinite(KL_div)) ).zeros();
  return(sum(KL_div));
}

// negative log-likelihood for log-prior variance, lV
double nll(const double lV, const vec& tau_no_V, const vec& nu, const vec& prior_weights) {
  vec tau = tau_no_V + exp(-lV);
  vec m = 0.5 * (-log(tau) + (pow(nu, 2) / tau));
  double m_max = max(m);
  vec w = exp(m - m_max);
  return((0.5 * lV) - m_max - log(sum(prior_weights % w)));
}

//double gr(const double lV, const vec& tau_no_V, const vec& nu, const vec& prior_weights) {
//    vec tau = tau_no_V + exp(-lV);
//    vec nu2_tau = pow(nu, 2) / tau;
//    vec m = 0.5 * (-log(tau) + nu2_tau);
//    double m_max = max(m);
//    vec w = exp(m - m_max);
//
//    double denom = dot(prior_weights, w) * exp(lV);
//    double numer = dot(prior_weights, (w / tau) % (1.0 + nu2_tau));
//
//    return(0.5 - 0.5 * (numer / denom));
//}

// [[Rcpp::export]]
List weighted_SER_cpp(XPtr<stumpsmatrix::StumpsMatrix> xp, arma::vec& Y, arma::vec& sigma2, Nullable<List> init = R_NilValue, double max_lV = 0, double lin_prior_prob = 0.5) {
  arma::vec s2 = arma::ones(Y.size());
  if (sigma2.size() == 1) {
    s2 = s2 * sigma2[0];
  } else {
    s2 = s2 % sigma2;
  }

  arma::vec inv_sigma2 = 1 / s2;
  double sum_inv_sigma2 = sum(inv_sigma2);
  arma::vec w = inv_sigma2 / sum_inv_sigma2;
  //double lV_init;
  arma::vec prior_weights;
  if (init.isNull() || !init.as().containsElementNamed("prior_weights") || !init.as().containsElementNamed("V")) {
    //lV_init = min(0.0, max_lV);
    double p = (double)xp.get()->ncol;
    double p_lin = (double)xp.get()->ncol_lin;
    double p_stumps = p - p_lin;
    prior_weights = arma::join_cols(arma::linspace(lin_prior_prob / p_lin, lin_prior_prob / p_lin, p_lin), arma::linspace((1 - lin_prior_prob) / p_stumps, (1 - lin_prior_prob) / p_stumps, p_stumps)) / (lin_prior_prob*(p_lin > 0) + (1 - lin_prior_prob)*(p_stumps > 0));
  } else {
    //lV_init = log((double)init.as()["V"]);
    prior_weights = as<arma::vec>(init.as()["prior_weights"]);
  }

  double Y_avg = arma::sum(Y % w);
  arma::vec Y_cent = Y - Y_avg;
  arma::vec X_avg = xp.get()->compute_Xty(w, arma::zeros(xp.get()->ncol));

  arma::vec tau_no_V = xp.get()->compute_X2ty(inv_sigma2, X_avg);
  arma::vec nu = xp.get()->compute_Xty(Y_cent / s2, X_avg);

  //double V = exp(optimize_lV(tau_no_V, nu, sigma2, prior_weights, max_lV, lV_init).at(0, 0));
  //Environment r_stats("package:stats");
  //Function r_optim = r_stats["optim"];
  //List optim_res = r_optim(Named("par") = lV_init, Named("fn") = InternalFunction(&nll), Named("tau_no_V") = tau_no_V, Named("nu") = nu,
  //                               Named("prior_weights") = prior_weights, Named("method") = "Brent", Named("lower") = -15, Named("upper") = max_lV);
  //arma::vec par = as<arma::vec>(optim_res["par"]);
  //double lV = optimlV::optimize_lV(lV_init, max_lV, Eigen::Map<Eigen::VectorXd>(tau_no_V.memptr(), tau_no_V.n_rows), 
  //                                                  Eigen::Map<Eigen::VectorXd>(nu.memptr(), nu.n_rows), 
  //                                                  Eigen::Map<Eigen::VectorXd>(prior_weights.memptr(), prior_weights.n_rows));
  //optimlV::NegLogLikR nll;
  //roptim::Roptim<optimlV::NegLogLikR> opt("L-BFGS-B");
  //opt.set_lower(-15.0);
  //opt.set_upper(max_lV);
  //opt.minimize(nll, init_lV);
  //double lV = opt.par()[0];
  Environment r_stats("package:stats");
  Function r_optimize = r_stats["optimize"];
  /*List optim_res = r_optim(Named("par") = lV_init, Named("fn") = InternalFunction(&nll), Named("gr") = InternalFunction(&gr), Named("tau_no_V") = tau_no_V, Named("nu") = nu,
                                Named("prior_weights") = prior_weights, Named("method") = "L-BFGS-B", Named("lower") = -15, Named("upper") = max_lV);*/
  List optim_res = r_optimize(Named("f") = InternalFunction(&nll),  Named("tau_no_V") = tau_no_V, Named("nu") = nu,
      Named("prior_weights") = prior_weights, Named("lower") = -15.0, Named("upper") = max_lV);
  arma::vec par = as<arma::vec>(optim_res["minimum"]);
  double V = std::exp(par[0]);
  //auto fn = [&, tau_no_V, nu, prior_weights](double lV) { return(nll(lV, tau_no_V, nu, prior_weights)); };
  //std::pair<double, double> opt = boost::math::tools::brent_find_minima(fn, -15.0, max_lV, 8);
  //double V = std::exp(opt.first);
  
  arma::vec tau = tau_no_V + (1 / V);

  arma::vec alpha = arma::log(prior_weights / arma::sqrt(tau)) + (0.5 * arma::pow(nu, 2) / tau);
  alpha = alpha - arma::max(alpha);
  alpha = arma::exp(alpha);
  alpha = alpha / arma::sum(alpha);

  arma::vec mu = nu / tau;
  
  arma::vec sigma2_post = 1 / tau;
  
  arma::vec beta_post_1 = alpha % mu;
  arma::vec beta_post_2 = alpha % (mu%mu + sigma2_post);

  arma::vec Xb_post = xp.get()->compute_Xb(beta_post_1, X_avg);

  arma::vec mu1 = Y_avg + Xb_post;
  arma::vec mu2 = std::pow(Y_avg, 2) + 2*Y_avg*Xb_post + xp.get()->compute_X2b(beta_post_2, X_avg);

  double kl_div = calc_KL_cpp(mu, alpha, sigma2_post, prior_weights, V);
  
  List res = List::create(Named("mu1") = mu1, Named("mu2") = mu2, Named("KL_div") = kl_div, Named("alpha") = alpha, Named("mu") = mu,
                          Named("sigma2_post") = sigma2_post, Named("V") = V, Named("X_avg") = X_avg, Named("Y_avg") = Y_avg, Named("prior_weights") = prior_weights);
  return(res);
}

// [[Rcpp::export]]
arma::vec predFnSusieStumps_cpp(XPtr<stumpsmatrix::StumpsMatrix> xp_new, List currentFit, int moment) {
  arma::vec beta_post_1 = as<arma::vec>(currentFit["alpha"]) % as<arma::vec>(currentFit["mu"]);
  arma::vec Xb_post = xp_new.get()->compute_Xb(beta_post_1, as<arma::vec>(currentFit["X_avg"]));
  if (moment == 1) {
    return(as<double>(currentFit["Y_avg"]) + Xb_post);
  } else if (moment == 2) {
    vec beta_post_2 = as<arma::vec>(currentFit["alpha"]) % (as<arma::vec>(currentFit["mu"])%as<arma::vec>(currentFit["mu"]) + as<arma::vec>(currentFit["sigma2_post"]));
    return(std::pow(as<double>(currentFit["Y_avg"]), 2) + 2*as<double>(currentFit["Y_avg"])*Xb_post + xp_new.get()->compute_X2b(beta_post_2, as<arma::vec>(currentFit["X_avg"])));
  } 
  stop("'moment' must be either 1 or 2");
}

