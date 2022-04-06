### Ridge Learner

#' @title Computes KL divergence for Ridge
#' @param mu_post is posterior mean of the coefficient vector
#' @param L_post is posterior precision of the coefficient vector
#' @param V is prior variance of the coefficient vector
calc_KL_Ridge_no_scale = function(mu_post, L_post, V) {
  p = length(mu_post)
  return(as.numeric(.5*(sum(1 / attr(L_post, 'evd')$values)/V + crossprod(mu_post)/V - p + (p*log(V)) + sum(log(attr(L_post, 'evd')$values)))))
}

#' @title negative log-likelihood for variance of ridge (call it V), on log-scale
#' @description b ~ N(0, VS) for fixed diagonal scaling matrix S
#' @description y = Xb + e
#' @description e ~ N(0, diag(sigma2)), L = diag(1/sigma2)
#' @param lV is ln(V)
#' @param XtY is X'Y
#' @param XtX is X'X (it must have an attribute called 'evd' that has its eigenvalue decomposition)
neg_log_lik_Ridge_logscale_no_scale = function(lV, XtY, XtX) {
  p = nrow(XtX)
  return(as.numeric(.5*p*lV + .5*sum(log(attr(XtX, 'evd')$values + exp(-lV))) - .5*sum(crossprod(attr(XtX, 'evd')$vectors, XtY)^2 / (attr(XtX, 'evd')$values + exp(-lV)))))
}


optimize_V_Ridge_no_scale = function(XtY, XtX, V = 1) {
  lV = optim(par = log(V), fn = neg_log_lik_Ridge_logscale_no_scale, XtY = XtY, XtX = XtX, method = 'Brent', lower = -15, upper = 35)$par
  V = exp(lV)
  return(V)
}


#' @title Ridge fit function
#' @param X is the design matrix. Scaling is not handled by the function, it must be dealt with by the user beforehand
#' @importFrom Matrix Diagonal
#' @importFrom Matrix crossprod
#' @importFrom Matrix tcrossprod
fit_ridge_no_scale = function(X, Y, sigma2, currentFit = list(V = NULL)) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(Y))
  }

  inv_sigma2 = 1 / sigma2
  sum_inv_sigma2 = sum(inv_sigma2)
  w = inv_sigma2 / sum_inv_sigma2

  # weighted center to deal with intercept
  Y_avg = sum(Y * w)
  Y_cent = Y - Y_avg
  X_avg = crossprod(X, w) # vector of weighted avg of columns of X
  X_cent = sweep(X, MARGIN = 2, STATS = X_avg, FUN = '-')

  # scale rows and response by 1/sqrt(sigma2) so we are in homoskedastic case
  Y_tilde = Y_cent / sqrt(sigma2)
  X_tilde = sweep(X_cent, MARGIN = 1, STATS = sqrt(sigma2), FUN = '/')

  XtY = as.numeric(crossprod(X_tilde, Y_tilde))
  XtX = crossprod(X_tilde)
  attr(XtX, 'evd') = eigen(XtX, symmetric = TRUE)

  if (is.null(currentFit$V)) {
    V = 1
  } else {
    V = currentFit$V
  }

  V = optimize_V_Ridge(XtY, XtX, V)

  L_post = XtX
  attr(L_post, 'evd')$values = attr(L_post, 'evd')$values + 1/V

  mu_post = tcrossprod(attr(L_post, 'evd')$vectors %*% Diagonal(x = 1 / sqrt(attr(L_post, 'evd')$values))) %*% XtY
  # XL_post_invXt = tcrossprod(X_tilde %*% attr(L_post, 'evd')$vectors %*% Diagonal(x = 1 / sqrt(attr(L_post, 'evd')$values)))

  Xb_post = X %*% mu_post
  X_avg_b_post = sum(X_avg * mu_post)

  intercept = as.numeric(Y_avg - X_avg_b_post)

  mu1 = as.numeric(Xb_post + intercept)
  # mu2 = Y_avg^2 + 2*Y_avg*(Xb_post - X_avg_b_post) + Matrix::diag(XL_post_invXt) + mu1^2
  mu2 = Y_avg^2 + 2*Y_avg*(Xb_post - X_avg_b_post) + apply(X_cent, MARGIN = 1, function(x) as.numeric(crossprod(Diagonal(x = 1 / sqrt(attr(L_post, 'evd')$values)) %*% crossprod(attr(L_post, 'evd')$vectors, x)) + sum(x * mu_post)^2))

  KL_div = calc_KL_Ridge(mu_post, L_post, V)

  return(list(mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), KL_div = KL_div, V = V, mu_post = as.numeric(mu_post), L_post_evd = attr(L_post, 'evd'), X_avg = X_avg, Y_avg = Y_avg, intercept = intercept))
}

fitFnRidge_no_scale = function(X, Y, sigma2, currentFit) {
  # supplied X here is the correlation  (called 'R' in the function above)
  return(fit_ridge_no_scale(X, Y, sigma2, currentFit))
}

# NOTE: this assumes that all obs in X_test are independent from obs used to fit
predFnRidge_no_scale = function(X_test, currentFit, moment = c(1, 2)) {
  Xb_post = as.numeric(X_test %*% currentFit$mu_post)
  if (moment == 1) {
    return(Xb_post + currentFit$intercept)
  } else if (moment == 2) {
    # return(Xb_post^2 + Matrix::diag(tcrossprod(X_test %*% L_post_evd$vectors %*% Diagonal(x = 1 / sqrt(L_post_evd$values)))))
    X_test_cent = sweep(X_test, MARGIN = 2, STATS = currentFit$X_avg, FUN = '-')
    return(currentFit$Y_avg^2 + 2*currentFit$Y_avg*(Xb_post - sum(currentFit$X_avg * currentFit$mu_post)) + apply(X_test_cent, MARGIN = 1, function(x) as.numeric(crossprod(Diagonal(x = 1 / sqrt(currentFit$L_post_evd$values)) %*% crossprod(currentFit$L_post_evd$vectors, x)) + sum(x * currentFit$mu_post)^2)))
  } else {
    stop("`moment` must be either 1 or 2")
  }
}


constCheckFnRidge_no_scale = function(currentFit) {
  return(currentFit$V < 1e-3)
}

#' Create a Ridge learner object
#'
#' Creates a Ridge learner object to be used in \code{\link{veb_boost}}
#'
#' @details A Ridge learner \eqn{\beta} has a prior distribution \deqn{\beta ~ N(0, \sigma_{\beta}^2 S^2)} for a fixed and given \eqn{S}.
#'
#' @param X is our design matrix used in training. Scaling is NOT handled by this function, and must be done by the user beforehand
#'
#' @param X_test is our design matrix used in testing Scaling is NOT handled by this function, and must be done by the user beforehand
#'
#' @param growMode is a string for if the learner should be grown (or not)
#' If \code{"+*"}, we grow mu_0 -> (mu_0 * mu_2) + mu_1
#' If \code{"+"}, we grow mu_0 -> (mu_0 + mu_1)
#' If \code{"*"}, we grow mu_0 -> (mu_0 * mu_1) (NOTE: Not recommended if we start with \code{k = 1})
#' If \code{"NA"}, we do not grow this learner
#'
#' @param changeToConstant is a logical for if constant fits should be changed to be constant
#'
#' @export
#'
makeRidgeLearner_no_scale = function(X, X_test = NULL, growMode = c("NA", "+*", "+", "*"), changeToConstant = FALSE) {
  growMode = match.arg(growMode)
  if (!(changeToConstant %in% c(TRUE, FALSE))) {
    stop("'changeToConstant' must be either TRUE or FALSE")
  }
  if (!is_valid_matrix(X)) {
    stop("'X' must be a numeric matrix")
  }
  if (!is.null(X_test) && !is_valid_matrix(X_test)) {
    stop("'X_test' must be a numeric matrix")
  }

  ridgeLearner = list(
    fitFunction = fitFnRidge_no_scale,
    predFunction = predFnRidge_no_scale,
    constCheckFunction = constCheckFnRidge_no_scale,
    currentFit = NULL,
    X = X,
    X_test = X_test,
    growMode = growMode,
    changeToConstant = changeToConstant
  )
  return(ridgeLearner)
}
