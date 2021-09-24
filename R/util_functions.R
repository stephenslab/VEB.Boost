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
#' @keywords internal
ilogit = function(x) { # inverse logit function
  1 / (1 + exp(-x))
}

# numerically stable log(1 + exp(x))
# see https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf formula 10
#' @keywords internal
log1pexp = function(x) {
  suppressWarnings(rowSums(cbind((exp(x) * (x <= -37)), 
                                 (log1p(exp(x)) * (x > -37) * (x <= 18)), 
                                 ((x + exp(-x)) * (x > 18) * (x <= 33.3)), 
                                 (x * (x > 33.3))), na.rm = T))
}

# numerically stable log(log(1 + exp(x)))
# see https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf formula 10
#' @keywords internal
loglog1pexp = function(x) {
  suppressWarnings(rowSums(cbind((x * (x <= -37)), 
                                 (log(log1p(exp(x))) * (x > -37) * (x <= 18)), 
                                 (log(x + exp(-x)) * (x > 18) * (x <= 33.3)), 
                                 (log(x) * (x > 33.3))), na.rm = T))
}
