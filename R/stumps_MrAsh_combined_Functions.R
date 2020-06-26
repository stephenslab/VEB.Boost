calc_ELBO = function(currentFit, Y, sigma2) { # calc ELBO up to constant term of -.5*sum(log(2*pi*sigma2))
  return(-.5 * (-2*sum((Y * currentFit$mu1 / sigma2)) + sum(currentFit$mu2 / sigma2)) - currentFit$KL_div)
}

fitFnStumpsMrAsh = function(X, Y, sigma2, init) {
  fit_Stumps = fitFnSusieStumps(X, Y, sigma2, init$fit_Stumps)
  fit_MrAsh = fitFn.mr.ash(X, Y, sigma2, init$fit_MrAsh)
  
  # p = ncol(X[[1]])
  # prior_weights = c((p - 1) / p, 1 / p)
  # prior_weights = c(.5, .5)
  p = get_ncol(X[-1])
  prior_weights = c(p / (p + 1), 1 / (p + 1))
  ELBOs = c(calc_ELBO(fit_Stumps, Y, sigma2), calc_ELBO(fit_MrAsh, Y, sigma2))

  post_weights = log(prior_weights) + ELBOs
  post_weights = post_weights - max(post_weights)
  post_weights = exp(post_weights)
  post_weights = post_weights / sum(post_weights)
  
  mu1 = (post_weights[1] * fit_Stumps$mu1) + (post_weights[2] * fit_MrAsh$mu1)
  mu2 = (post_weights[1] * fit_Stumps$mu2) + (post_weights[2] * fit_MrAsh$mu2)
  KL_div = post_weights * (log(post_weights) - log(prior_weights) + c(fit_Stumps$KL_div, fit_MrAsh$KL_div))
  KL_div[post_weights == 0] = 0
  KL_div = sum(KL_div)
  
  return(list(mu1 = mu1, mu2 = mu2, KL_div = KL_div, prior_weights = prior_weights, post_weights = post_weights, fit_Stumps = fit_Stumps, fit_MrAsh = fit_MrAsh))
}

predFnStumpsMrAsh = function(X_new, currentFit, moment = c(1, 2)) {
  pred_Stumps = predFnSusieStumps(X_new[-1], currentFit$fit_Stumps, moment)
  pred_MrAsh = predFn.mr.ash(X_new, currentFit$fit_MrAsh, moment)
  
  return((currentFit$post_weights[1] * pred_Stumps) + (currentFit$post_weights[2] * pred_MrAsh))
}

constCheckFnStumpsMrAsh = function(currentFit) {
  V_Stumps = currentFit$fit_Stumps$V
  V_MrAsh = sum(currentFit$fit_MrAsh$prior_pi * currentFit$fit_MrAsh$prior_var)
  V = (currentFit$prior_weights[1] * V_Stumps) + (currentFit$prior_weights[2] * V_MrAsh)
  return(V < 1e-3)
  # return(max(currentFit$mu2 - currentFit$mu1^2) < 1e-6)
}
