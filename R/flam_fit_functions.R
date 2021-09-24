### FLAM-related functions

# weighted FLAM function, ONLY stumps
# X is a list, only with stumps for variables
weighted_flam = function(X, Y, sigma2, init = list(flam_learners = NULL)) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(Y))
  }
  
  # can parallelize easily with mclapply
  if (is.null(init$flam_learners)) {
    flam_learners = lapply(X, function(x) {
      learner_x = veb_boost(list(x), Y, fitFunctions = fitFnSusieStumps, predFunctions = predFnSusieStumps, constCheckFunctions = constCheckFnSusieStumps,
                            k = 10, growTree = FALSE, sigma2 = sigma2, family = "gaussian", verbose = FALSE, changeToConstant = FALSE)
      return(learner_x)
    })
  } else {
    flam_learners = lapply(init$flam_learners, function(x) {
      x$Y = Y
      x$sigma2 = sigma2
      x$convergeFit(tol = length(Y) / 10000, update_sigma2 = FALSE, update_ELBO_progress = FALSE, verbose = FALSE, maxit = 10)
      return(x)
    })
  }

  p = length(flam_learners)
  prior_weights = rep(1 / p, p)
  ELBOs = sapply(flam_learners, function(x) x$ELBO)
  
  post_weights = log(prior_weights) + ELBOs
  post_weights = post_weights - max(post_weights)
  post_weights = exp(post_weights)
  post_weights = post_weights / sum(post_weights)
  
  mu1 = rowSums(sweep(do.call(cbind, lapply(flam_learners, function(x) x$mu1)), 2, post_weights, '*'))
  mu2 = rowSums(sweep(do.call(cbind, lapply(flam_learners, function(x) x$mu2)), 2, post_weights, '*'))
  KL_div = post_weights * (log(post_weights) - log(prior_weights) + sapply(flam_learners, function(x) x$KL_div))
  KL_div[post_weights == 0] = 0
  KL_div = sum(KL_div)

  return(list(mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), KL_div = KL_div, flam_learners = flam_learners, post_weights = post_weights))
}

fitFnFLAM = function(X, Y, sigma2, init) {
  return(weighted_flam(X, Y, sigma2, init))
}

predFnFLAM = function(X_new, currentFit, moment = c(1, 2)) {
  lapply(currentFit$flam_learners, function(x) x$predict(X_new, moment))
  if (moment == 1) {
    return(rowSums(sweep(do.call(cbind, lapply(currentFit$flam_learners, function(x) x$pred_mu1)), 2, currentFit$post_weights, '*')))
  } else if (moment == 2) {
    return(rowSums(sweep(do.call(cbind, lapply(currentFit$flam_learners, function(x) x$pred_mu2)), 2, currentFit$post_weights, '*')))
  } else {
    stop("`moment` must be either 1 or 2")
  }
}


constCheckFnFLAM = function(currentFit) {
  V_maxs = sapply(currentFit$flam_learners, function(x) max(sapply(x$leaves, function(x) x$currentFit$V), na.rm = TRUE))
  return(max(V_maxs) < 1e-3)
}
