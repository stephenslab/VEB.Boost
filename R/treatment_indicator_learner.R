#' @keywords internal
fitFnTrt = function(X, Y, sigma2, currentFit) { # treatment indicator learner fit function
  KL_div = 0

  mu1 = currentFit$mu1
  mu2 = currentFit$mu2
  return(list(mu1 = mu1, mu2 = mu2, KL_div = KL_div))
}
#' @keywords internal
predFnTrt = function(X_test, currentFit, moment = c(1, 2)) { # treatment indicator learner prediction function
  if (moment == 1) {
    return(currentFit$mu1)
  } else if (moment == 2) {
    return(currentFit$mu2)
  } else {
    stop("`moment` must be either 1 or 2")
  }
}
#' @keywords internal
constCheckFnTrt = function(currentFit) { # treatment indicator learner constant check function
  return(TRUE)
}

# treatment indicator learner
#' @keywords internal
makeTrtLearner = function(z) {
  trtLearner = list(
    fitFunction = fitFnTrt,
    predFunction = predFnTrt,
    constCheckFunction = constCheckFnTrt,
    currentFit = list(mu1 = z, mu2 = z, KL_div = 0),
    X = NULL,
    X_test = NULL,
    growMode = "NA",
    changeToConstant = FALSE
  )

  return(trtLearner)
}
