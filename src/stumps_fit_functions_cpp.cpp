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
  vec KL_div = alpha % (log((alpha / prior_weights) % sqrt(V / sigma2)) - .5 + ((sigma2 + (mu % mu)) / (2 * V)));
  KL_div.elem( find(alpha == 0.0) ).zeros();
  // KL_div.elem( find(approx_equal(alpha, 0.0, "absdiff", 1e-12)) ) = 0.0;
  // KL_div.elem( find(find_nonfinite(KL_div)) ).zeros();
  return(sum(KL_div));
    /*
  double res = 0;
  double cumsum_alpha = 0;
  //for (size_t i = 0; i < KL_div.size(); i++) {
  for (size_t i = 0; i < alpha.size(); i++) {
      if (alpha[i] != 0.0) {
          //res += KL_div[i];
          res += alpha[i] * (log((alpha[i] / prior_weights[i]) * sqrt(V / sigma2[i])) - .5 + ((sigma2[i] + (mu[i] * mu[i])) / (2 * V)));
          cumsum_alpha += alpha[i];
          if (cumsum_alpha == 1.0) {
              break;
          }
      }
  }
  return(res);
  */
}

// negative log-likelihood for log-prior variance, lV
double nll(const double lV, const vec& tau_no_V, const vec& nu2, const vec& prior_weights) {
  vec tau = tau_no_V + exp(-lV);
  vec m = 0.5 * (-arma::log(tau) + (nu2 / tau));
  double m_max = arma::max(m);
  vec w = arma::exp(m - m_max);
  return((0.5 * lV) - m_max - log(arma::sum(prior_weights % w)));
}

/*
class nllFunctorClass : public brent::func_base
{
private:
    const vec& tau_no_V, nu2, prior_weights;
public:
    nllFunctorClass(const vec& tau_no_V_, const vec& nu2_, const vec& prior_weights_) : tau_no_V(tau_no_V_), nu2(nu2_), prior_weights(prior_weights_) {}
    double operator() (double lV) {
        vec tau = tau_no_V + exp(-lV);
        vec m = 0.5 * (-arma::log(tau) + (nu2 / tau));
        double m_max = arma::max(m);
        vec w = arma::exp(m - m_max);
        return((0.5 * lV) - m_max - log(arma::sum(prior_weights % w)));
    }
};
*/
/*
double gr(const double lV, const vec& tau_no_V, const vec& nu2, const vec& prior_weights) {
    vec tau = tau_no_V + exp(-lV);
    vec nu2_tau = nu2 / tau;
    vec m = 0.5 * (-log(tau) + nu2_tau);
    double m_max = max(m);
    vec w = exp(m - m_max);

    double denom = sum(prior_weights % w) * exp(lV);
    double numer = sum(prior_weights % (w / tau) % (1.0 + nu2_tau));

    return(0.5 - 0.5 * (numer / denom));
}
*/

// [[Rcpp::export]]
List weighted_SER_cpp(XPtr<stumpsmatrix::StumpsMatrix> xp, arma::vec& Y, arma::vec& sigma2, Nullable<List> init = R_NilValue, double max_lV = 0.0, double lin_prior_prob = 0.5, bool use_optim = true) {
  arma::vec s2 = arma::ones(Y.size());
  if (sigma2.size() == 1) {
    s2 *= sigma2[0];
  } else {
    s2 %= sigma2;
  }

  arma::vec inv_sigma2 = 1 / s2;
  inv_sigma2.replace(arma::datum::nan, 0.0);
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

   arma::mat X_avg_tau_no_V_nu = xp.get()->compute_X_avg_tau_no_V_nu(w, inv_sigma2, Y_cent);
   arma::vec nu2 = X_avg_tau_no_V_nu.col(2) % X_avg_tau_no_V_nu.col(2);
   //arma::vec X_avg = X_avg_tau_no_V_nu.col(0);
   //arma::vec tau_no_V = X_avg_tau_no_V_nu.col(1);
   //arma::vec nu = X_avg_tau_no_V_nu.col(2);
  
  /*
  arma::vec X_avg = xp.get()->compute_Xty(w, arma::zeros(xp.get()->ncol));

  arma::vec tau_no_V = xp.get()->compute_X2ty(inv_sigma2, X_avg);
  arma::vec nu = xp.get()->compute_Xty(Y_cent % inv_sigma2, X_avg);
  */

  //double V = exp(optimize_lV(tau_no_V, nu, sigma2, prior_weights, max_lV, lV_init).at(0, 0));
  //Environment r_stats("package:stats");
  //Function r_optim = r_stats["optim"];
  //List optim_res = r_optim(Named("par") = lV_init, Named("fn") = InternalFunction(&nll), Named("tau_no_V") = tau_no_V, Named("nu2") = nu2,
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
  double V;
  double lV_init = std::min(0.0, max_lV);
  if (init.as().containsElementNamed("V")) {
      lV_init = log((double)init.as()["V"]);
  }
  if (use_optim || init.isNull() || !init.as().containsElementNamed("V")) { // if using Brent method to optimize for V....
      Environment r_stats("package:stats");
      Function r_optimize = r_stats["optimize"];
      List optim_res = r_optimize(Named("f") = InternalFunction(&nll), Named("tau_no_V") = X_avg_tau_no_V_nu.col(1), Named("nu2") = nu2,
          Named("prior_weights") = prior_weights, Named("lower") = -15.0, Named("upper") = max_lV);
      /*
      Function r_optim = r_stats["optim"];
      List optim_res = r_optim(Named("par") = lV_init, Named("fn") = InternalFunction(&nll), Named("gr") = InternalFunction(&gr), Named("tau_no_V") = X_avg_tau_no_V_nu.col(1), Named("nu2") = nu2,
                                    Named("prior_weights") = prior_weights, Named("method") = "L-BFGS-B", Named("lower") = -15, Named("upper") = max_lV);
          */
      arma::vec par = as<arma::vec>(optim_res["minimum"]);
      V = std::exp(par[0]);
      /*
      auto fn = [&, tau_no_V, nu, prior_weights](double lV) { return(nll(lV, tau_no_V, nu, prior_weights)); };
      std::pair<double, double> opt = boost::math::tools::brent_find_minima(fn, -15.0, max_lV, 8);
      double V = std::exp(opt.first);
      */
      /*
      nllFunctorClass f(X_avg_tau_no_V_nu.col(1), X_avg_tau_no_V_nu.col(2) % X_avg_tau_no_V_nu.col(2), prior_weights);
      double f_star = brent::local_min(-15.0, max_lV, 0.0001220703, f, V); // 0.0001220703 = .Machine$double.eps^0.25 = std::pow(DBL_EPSILON, 0.25)
      V = std::exp(V);      
      */
  } else { // otherwise, using EM method
      V = (double)init.as()["V"];
  }
  
  arma::vec tau = X_avg_tau_no_V_nu.col(1) + (1 / V);

  arma::vec alpha = arma::log(prior_weights / arma::sqrt(tau)) + (0.5 * nu2 / tau);
  alpha -= arma::max(alpha);
  alpha = arma::exp(alpha); //alpha.transform([](double val) { return(std::exp(val)); });
  alpha /= arma::sum(alpha);

  arma::vec mu = X_avg_tau_no_V_nu.col(2) / tau;
  
  arma::vec sigma2_post = 1 / tau;

  if (!use_optim) { // if using EM method for V, update now
      V = arma::sum(alpha % (sigma2_post + (mu % mu)));
      V = std::max(std::exp(-15.0), std::min(V, std::exp(max_lV)));
  }
  
  arma::vec beta_post_1 = alpha % mu;
  arma::vec beta_post_2 = alpha % ((mu % mu) + sigma2_post);

  arma::vec Xb_post = xp.get()->compute_Xb(beta_post_1, X_avg_tau_no_V_nu.col(0));

  arma::vec mu1 = Y_avg + Xb_post;
  arma::vec mu2 = std::pow(Y_avg, 2) + 2*Y_avg*Xb_post + xp.get()->compute_X2b(beta_post_2, X_avg_tau_no_V_nu.col(0));

  double kl_div = calc_KL_cpp(mu, alpha, sigma2_post, prior_weights, V);
  
  List res = List::create(Named("mu1") = mu1, Named("mu2") = mu2, Named("KL_div") = kl_div, Named("alpha") = alpha, Named("mu") = mu,
                          Named("sigma2_post") = sigma2_post, Named("V") = V, Named("X_avg") = X_avg_tau_no_V_nu.col(0), Named("Y_avg") = Y_avg, Named("prior_weights") = prior_weights);
  return(res);
}

// [[Rcpp::export]]
arma::vec predFnSusieStumps_cpp(XPtr<stumpsmatrix::StumpsMatrix> X_test, List currentFit, int moment) {
  arma::vec beta_post_1 = as<arma::vec>(currentFit["alpha"]) % as<arma::vec>(currentFit["mu"]);
  arma::vec Xb_post = X_test.get()->compute_Xb(beta_post_1, as<arma::vec>(currentFit["X_avg"]));
  if (moment == 1) {
    return(as<double>(currentFit["Y_avg"]) + Xb_post);
  } else if (moment == 2) {
    vec beta_post_2 = as<arma::vec>(currentFit["alpha"]) % (as<arma::vec>(currentFit["mu"])%as<arma::vec>(currentFit["mu"]) + as<arma::vec>(currentFit["sigma2_post"]));
    return(std::pow(as<double>(currentFit["Y_avg"]), 2) + 2*as<double>(currentFit["Y_avg"])*Xb_post + X_test.get()->compute_X2b(beta_post_2, as<arma::vec>(currentFit["X_avg"])));
  } 
  stop("'moment' must be either 1 or 2");
}


// [[Rcpp::export]]
arma::vec getAlphaByVar(XPtr<stumpsmatrix::StumpsMatrix> xp, List currentFit) { // combine all probability for same variable (i.e. combine linear and stumps probabilities)
    arma::uvec ncol_vec = xp.get()->get_ncol_vec();
    NumericVector alpha = currentFit["alpha"];
    arma::uvec which_lin = arma::find(xp.get()->include_linear);
    arma::uvec which_stumps = arma::find(xp.get()->include_stumps);
    arma::vec alpha_comb(xp.get()->include_linear.size(), arma::fill::zeros);
    unsigned int offset = 0; 
    unsigned int j = 0;
    bool any_lin = xp.get()->ncol_lin > 0;

    if (any_lin) { // add linear probs, if needed
        for (size_t k = 0; k < which_lin.size(); k++) {
            alpha_comb[which_lin[k]] += alpha[offset];
            offset++;
        }
        //alpha_comb.elem(which_lin) += alpha[Range(offset, offset + ncol_vec[j] - 1)];
        //offset += ncol_vec[j];
        j++;
    }
    
    while (j < ncol_vec.size()) { // add stumps probs, if needed
        alpha_comb[which_stumps[j - any_lin]] += sum(alpha[Range(offset, offset + ncol_vec[j] - 1)]);
        offset += ncol_vec[j];
        j++;
    }

    return alpha_comb;
}
