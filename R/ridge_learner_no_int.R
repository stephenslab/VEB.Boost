### Ridge Learner

#' @title Computes KL divergence for Ridge
#' @param mu_post is posterior mean of the coefficient vector
#' @param L_post is posterior precision of the coefficient vector
#' @param V is prior variance of the coefficient vector
calc_KL_Ridge_no_int = function(mu_post, L_post, V) {
  p = length(mu_post)
  return(as.numeric(.5*(sum(1 / attr(L_post, 'evd')$values)/V + crossprod(mu_post)/V - p + (p*log(V)) + sum(log(attr(L_post, 'evd')$values)))))
}

#' @title negative log-likelihood for variance of ridge (call it V), on log-scale
#' @description b ~ N(0, VS) for fixed diagonal scaling matrix S
#' @description y = Xb + e
#' @description e ~ N(0, diag(sigma2)), L = diag(1/sigma2)
#' @param lV is ln(V)
#' @param XtLY is X'LY
#' @param XtLX is X'LX (it must have an attribute called 'evd' that has its eigenvalue decomposition)
neg_log_lik_Ridge_logscale_no_int = function(lV, XtLY, XtLX) {
  p = nrow(XtLX)
  return(as.numeric(.5*p*lV + .5*sum(log(attr(XtLX, 'evd')$values + exp(-lV))) - .5*sum(crossprod(attr(XtLX, 'evd')$vectors, XtLY)^2 / (attr(XtLX, 'evd')$values + exp(-lV)))))
}


optimize_V_Ridge_no_int = function(XtLY, XtLX, V = 1) {
  lV = optim(par = log(V), fn = neg_log_lik_Ridge_logscale_no_int, XtLY = XtLY, XtLX = XtLX, method = 'Brent', lower = -15, upper = 35)$par
  V = exp(lV)
  return(V)
}


#' @title Ridge fit function
#' @param X is the design matrix. Scaling is not handled by the function, it must be dealt with by the user beforehand
#' @importFrom Matrix Diagonal
#' @importFrom Matrix crossprod
#' @importFrom Matrix tcrossprod
fit_ridge_no_int = function(X, Y, sigma2, init = list(V = NULL)) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(Y))
  }
  L = Diagonal(x = 1 / sigma2)

  XtLY = as.numeric(crossprod(X, L %*% Y))
  XtLX = crossprod(sqrt(L) %*% X)
  attr(XtLX, 'evd') = eigen(XtLX, symmetric = TRUE)

  if (is.null(init$V)) {
    V = 1
  } else {
    V = init$V
  }

  V = optimize_V_Ridge_no_int(XtLY, XtLX, V)

  L_post = XtLX
  attr(L_post, 'evd')$values = attr(L_post, 'evd')$values + 1/V

  mu_post = tcrossprod(attr(L_post, 'evd')$vectors %*% Diagonal(x = 1 / sqrt(attr(L_post, 'evd')$values))) %*% XtLY
  XL_post_invXt = tcrossprod(X %*% attr(L_post, 'evd')$vectors %*% Diagonal(x = 1 / sqrt(attr(L_post, 'evd')$values)))

  mu1 = as.numeric(XL_post_invXt %*% L %*% Y)
  mu2 = mu1^2 + Matrix::diag(XL_post_invXt)

  KL_div = calc_KL_Ridge(mu_post, L_post, V)

  return(list(mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), KL_div = KL_div, V = V, mu_post = as.numeric(mu_post), L_post_evd = attr(L_post, 'evd')))
}

fitFnRidge_no_int = function(X, Y, sigma2, init) {
  # supplied X here is the correlation  (called 'R' in the function above)
  return(fit_ridge_no_int(X, Y, sigma2, init))
}

# NOTE: this assumes that all obs in X_new are independent from obs used to fit
predFnRidge_no_int = function(X_new, currentFit, moment = c(1, 2)) {
  Xb_post = as.numeric(X_new %*% currentFit$mu_post)
  if (moment == 1) {
    return(Xb_post)
  } else if (moment == 2) {
    return(Xb_post^2 + Matrix::diag(tcrossprod(X_new %*% L_post_evd$vectors %*% Diagonal(x = 1 / sqrt(L_post_evd$values)))))
  } else {
    stop("`moment` must be either 1 or 2")
  }
}


constCheckFnRidge_no_int = function(currentFit) {
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
makeRidgeLearner_no_int = function(X, X_test = NULL, growMode = c("NA", "+*", "+", "*"), changeToConstant = FALSE) {
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
    fitFunction = fitFnRidge_no_int,
    predFunction = predFnRidge_no_int,
    constCheckFunction = constCheckFnRidge_no_int,
    currentFit = NULL,
    X = X,
    X_test = X_test,
    growMode = growMode,
    changeToConstant = changeToConstant
  )
  return(ridgeLearner)
}
