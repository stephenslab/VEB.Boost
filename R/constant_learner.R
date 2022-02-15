#' @keywords internal
fitFnConstComp = function(X, Y, sigma2, init) { # constant fit function
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(Y))
  }
  intercept = weighted.mean(Y, 1/sigma2)
  KL_div = 0
  
  mu1 = intercept
  mu2 = intercept^2
  return(list(mu1 = mu1, mu2 = mu2, intercept = intercept, KL_div = KL_div))
}
#' @keywords internal
predFnConstComp = function(X_new, currentFit, moment = c(1, 2)) { # constant prediction function
  if (moment == 1) {
    return(currentFit$intercept)
  } else if (moment == 2) {
    return(currentFit$intercept^2)
  } else {
    stop("`moment` must be either 1 or 2")
  }
}
#' @keywords internal
constCheckFnConstComp = function(currentFit) { # constant constant check function
  return(TRUE)
}

# constant learner
#' @keywords internal
constLearner = list(
  fitFunction = fitFnConstComp,
  predFunction = predFnConstComp,
  constCheckFunction = constCheckFnConstComp,
  currentFit = NULL,
  X = NULL,
  X_test = NULL,
  growMode = NULL,
  changeToConstant = FALSE
)