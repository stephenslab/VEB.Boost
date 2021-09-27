### SuSiE stumps-related functions

calc_KL = function(mu, alpha, sigma2, prior_weights, V = 1) {
  p = length(mu)
  b_post = alpha * mu

  prior_var = rep(V, p)

  KL_div = alpha * (log(alpha) - log(prior_weights) + (log(prior_var) / 2) - (log(sigma2) / 2) - .5 + ((sigma2 + mu^2) / (2 * prior_var)))
  KL_div[alpha == 0] = 0
  return(sum(KL_div))
}

# over-write optimize_V code w/ weighted version
log_lik_SER = function(V, tau_no_V, nu, sigma2, prior_weights) {
  tau = tau_no_V + (1 / V)
  m = -(log(tau) / 2) + (nu^2 / (2 * tau))
  m_max = max(m)
  w = exp(m - m_max)
  -(log(V) / 2) + m_max + log(sum(prior_weights * w))
}

neg.loglik.logscale = function(lV, tau_no_V, nu, sigma2, prior_weights){
  -log_lik_SER(exp(lV), tau_no_V, nu, sigma2, prior_weights)
}

optimize_V = function(tau_no_V, nu, sigma2, prior_weights, lV = min(0, max_lV), max_lV = 0) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(nu))
  }
  lV = optim(par = lV, fn = neg.loglik.logscale, tau_no_V = tau_no_V, nu = nu, sigma2 = sigma2, prior_weights = prior_weights, method='Brent', lower = -15, upper = max_lV)$par
  V = exp(lV)
  return(V)
}



# weighted SER function, linear terms + stumps
# X is a list, first element corresponds to linear, others are stumps for variables
weighted_SER = function(X, Y, sigma2, init = list(V = NULL), max_lV = 0, lin_prior_prob = 0.5) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(Y))
  }

  inv_sigma2 = 1 / sigma2
  sum_inv_sigma2 = sum(inv_sigma2)
  w = inv_sigma2 / sum_inv_sigma2
  p = get_ncol(X)
  p_lin = 0
  if (is_valid_matrix(X[[1]])) { # if the first matrix is linear terms
    p_lin = get_ncol(X[[1]])
  }
  p_stumps = p - p_lin
  # prior_weights = c(rep(.5 / p_lin, p_lin), rep(.5 / p_stumps, p_stumps)) * ifelse(p_lin * p_stumps == 0, 2, 1)
  prior_weights = c(rep(lin_prior_prob / p_lin, p_lin), rep((1 - lin_prior_prob) / p_stumps, p_stumps)) / (lin_prior_prob*(p_lin != 0) + (1 - lin_prior_prob)*(p_stumps != 0))
  # prior_weights = rep(1 / p, p)

  Y_avg = sum(Y * w)
  Y_cent = Y - Y_avg
  X_avg = compute_Xty(X, w) # vector of weighted avg of columns of X

  tau_no_V = compute_X2ty(X, inv_sigma2, X_avg)
  nu = compute_Xty(X, Y_cent / sigma2) - (X_avg * sum(Y_cent / sigma2))

  # optim method, seems to be slower than EM method
  lV = log(ifelse(is.null(init$V), 1, init$V))
  V = optimize_V(tau_no_V, nu, sigma2, prior_weights, lV, max_lV)

  tau = tau_no_V + (1 / V)

  alpha = log(prior_weights) - (.5 * log(tau)) + (.5 * nu^2 / tau)
  alpha = alpha - max(alpha)
  alpha = exp(alpha)
  alpha = alpha / sum(alpha)

  mu = nu / tau

  sigma2_post = 1 / tau

  # iterative EM version, seems to be faster than optim method (but sometimes takes HOURS to converve.... probably not a great idea)
  # V = ifelse(is.null(init$V), 1, init$V)
  # V_old = Inf
  # while(abs(V - V_old) > 1e-10) {
  #   V_old = V
  #   tau = tau_no_V + (1 / V)
  #
  #   alpha = log(prior_weights) - (.5 * log(tau)) + (.5 * nu^2 / tau)
  #   alpha = alpha - max(alpha)
  #   alpha = exp(alpha)
  #   alpha = alpha / sum(alpha)
  #
  #   mu = nu / tau
  #
  #   sigma2_post = 1 / tau
  #   V = sum(alpha * (sigma2_post + mu^2))
  # }

  # single EM update
  # V = ifelse(is.null(init$V), 1, init$V)
  # tau = tau_no_V + (1 / V)
  # alpha = log(prior_weights) - (.5 * log(tau)) + (.5 * nu^2 / tau)
  # alpha = alpha - max(alpha)
  # alpha = exp(alpha)
  # alpha = alpha / sum(alpha)
  #
  # mu = nu / tau
  #
  # sigma2_post = 1 / tau
  # V = sum(alpha * (sigma2_post + mu^2))

  beta_post_1 = alpha * mu
  beta_post_2 = alpha * (sigma2_post + mu^2)

  Xb_post = compute_Xb(X, beta_post_1)
  X_avg_b_post = sum(X_avg * beta_post_1)

  intercept = as.numeric(Y_avg - X_avg_b_post)

  # mu1 = E[int + Xb] = E[Y_avg - X_avg'b + Xb]
  mu1 = intercept + Xb_post
  # mu2 = E[(int + Xb)^2] = E[(Y_avg - X_avg'b + Xb)^2]
  mu2 = Y_avg^2 + 2*Y_avg*(Xb_post - X_avg_b_post) + compute_X2b(X, beta_post_2, X_avg)

  KL_div = calc_KL(mu, alpha, sigma2_post, prior_weights, V)

  return(list(mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), KL_div = KL_div, alpha = alpha, mu = mu, sigma2_post = sigma2_post, intercept = intercept, V = V, X_avg = X_avg, Y_avg = Y_avg))
}

fitFnSusieStumps = function(X, Y, sigma2, init) {
  return(weighted_SER(X, Y, sigma2, init, 0, 0.5))
}

predFnSusieStumps = function(X_new, currentFit, moment = c(1, 2)) {
  beta_post_1 = currentFit$alpha * currentFit$mu
  if (moment == 1) {
    return(currentFit$intercept + compute_Xb(X_new, beta_post_1))
  } else if (moment == 2) {
    beta_post_2 = currentFit$alpha * (currentFit$sigma2_post + currentFit$mu^2)
    return(currentFit$Y_avg^2 + 2*currentFit$Y_avg*(compute_Xb(X_new, beta_post_1) - sum(currentFit$X_avg * beta_post_1)) + compute_X2b(X_new, beta_post_2, currentFit$X_avg))
  } else {
    stop("`moment` must be either 1 or 2")
  }
}


constCheckFnSusieStumps = function(currentFit) {
  return(currentFit$V < 1e-3)
}
