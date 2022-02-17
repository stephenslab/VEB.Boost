### Random Effect (RE) related functions

#' @title Computes KL divergence for RE
#' @param mu_post is posterior mean of RE
#' @param L_post is posterior precision of RE (it must have an attribute called 'evd' that has its eigenvalue decomposition)
#' @param V is prior variance of RE
calc_KL_RE = function(mu_post, L_post, V) {
  k = length(mu_post)
  # chol_L_post_inv = backsolve(chol_L_post, diag(1, nrow(chol_L_post), ncol(chol_L_post)), transpose = TRUE)
  #
  # return(as.numeric(.5*(sum(chol_L_post_inv[upper.tri(chol_L_post_inv, diag = TRUE)]^2)/V + crossprod(mu_post)/V - k + (k*log(V)) + 2*sum(log(diag(chol_L_post))))))
  return(as.numeric(.5*(sum(1 / attr(L_post, 'evd')$values)/V + crossprod(mu_post)/V - k + (k*log(V)) + sum(log(attr(L_post, 'evd')$values)))))
}

#' @title negative log-likelihood for variance of random effect (call it V), on log-scale
#' @description a ~ N(0, VR) for fixed correlation matrix R = USU' (of rank k <= n)
#' @description z ~ N(0, VI_k)
#' @description a = US^(1/2)z
#' @description y = a + e
#' @description e ~ N(0, diag(sigma2)), L = diag(1/sigma2)
#' @param lV is ln(V)
#' @param SUtLY is S^(1/2)U'LY
#' @param SUtLUS is S^(1/2)U'LUS^(1/2) (it must have an attribute called 'evd' that has its eigenvalue decomposition)
neg_log_lik_RE_logscale = function(lV, SUtLY, SUtLUS) {
  k = nrow(SUtLUS)
  # M = as.matrix(SUtLUS)
  # Matrix::diag(M) = Matrix::diag(M) + rep(exp(-lV), k)
  # chol_M = chol(M)
  # return(as.numeric(.5*k*lV + sum(log(diag(chol_M))) - .5*crossprod(backsolve(chol_M, SUtLY, transpose = TRUE))))
  return(as.numeric(.5*k*lV + .5*sum(log(attr(SUtLUS, 'evd')$values + exp(-lV))) - .5*sum((attr(SUtLUS, 'evd')$vectors %*% SUtLY)^2 / (attr(SUtLUS, 'evd')$values + exp(-lV)))))
}


optimize_V_RE = function(SUtLY, SUtLUS, V = 1) {
  lV = optim(par = log(V), fn = neg_log_lik_RE_logscale, SUtLY = SUtLY, SUtLUS = SUtLUS, method='Brent', lower = -15, upper = 35)$par
  V = exp(lV)
  return(V)
}


#' @title RE fit function
#' @param R is the fixed correlation matrix (needs attribute 'svd' which has the SVD saved (with elements 'd' and 'u'))
#' @importFrom Matrix Diagonal
#' @importFrom Matrix crossprod
#' @importFrom Matrix tcrossprod
fit_random_effect = function(R, Y, sigma2, init = list(V = NULL)) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(Y))
  }
  U = attr(R, 'svd')$u
  S = Diagonal(x = sqrt(attr(R, 'svd')$d))
  L = Diagonal(x = 1 / sigma2)

  SUtLY = as.numeric(S %*% crossprod(U, L %*% Y))
  SUtLUS = S %*% crossprod(U, crossprod(L, U)) %*% S
  if ((nrow(U) == ncol(U)) && all(U == Matrix::Diagonal(nrow(U))) && all(attr(R, 'svd')$d == 1)) { # if U and S are both the identity....
    attr(SUtLUS, 'evd') = list(vectors = U, values = 1 / sigma2)
  } else {
    attr(SUtLUS, 'evd') = eigen(SUtLUS, symmetric = TRUE)
  }

  if (is.null(init$V)) {
    V = 1
  } else {
    V = init$V
  }

  V = optimize_V_RE(SUtLY, SUtLUS, V)

  L_post = SUtLUS
  attr(L_post, 'evd')$values = attr(L_post, 'evd')$values + 1/V
  # Matrix::diag(L_post) = Matrix::diag(L_post) + rep(1 / V, nrow(S))
  # chol_L_post = chol(L_post)
  # mu_post = as.numeric(solve(L_post, SUtLY))
  # mu_post = as.numeric(backsolve(chol_L_post, backsolve(chol_L_post, SUtLY, transpose = TRUE)))
  # USL_post_invSUt = crossprod(backsolve(chol_L_post, tcrossprod(S, U), transpose = TRUE))
  mu_post = crossprod(attr(L_post, 'evd')$vectors, (attr(L_post, 'evd')$vectors %*% SUtLY) / attr(L_post, 'evd')$values)
  USL_post_invSUt = crossprod(sweep(attr(L_post, 'evd')$vectors %*% tcrossprod(S, U), 1, sqrt(attr(L_post, 'evd')$values), '/'))

  mu1 = as.numeric(USL_post_invSUt %*% L %*% Y)
  # chol_L_post_inv =  backsolve(chol_L_post, Matrix::diag(rep(1, nrow(chol_L_post))))
  mu2 = mu1^2 + Matrix::diag(USL_post_invSUt)

  KL_div = calc_KL_RE(mu_post, L_post, V)

  return(list(mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), KL_div = KL_div, V = V))
}

fitFnRandomEffect = function(X, Y, sigma2, init) {
  # supplied X here is the correlation  (called 'R' in the function above)
  return(fit_random_effect(X, Y, sigma2, init))
}

# NOTE: this assumes that all obs in X_new are independent from obs used to fit
predFnRandomEffect = function(X_new, currentFit, moment = c(1, 2)) {
  beta_post_1 = currentFit$alpha * currentFit$mu
  if (moment == 1) {
    return(0)
  } else if (moment == 2) {
    return(V)
  } else {
    stop("`moment` must be either 1 or 2")
  }
}


constCheckFnRandomEffect = function(currentFit) {
  return(currentFit$V < 1e-3)
}

# constant learner
#' @export
randomEffectLearner = list(
  fitFunction = fitFnRandomEffect,
  predFunction = predFnRandomEffect,
  constCheckFunction = constCheckFnRandomEffect,
  currentFit = NULL,
  X = NULL,
  X_test = NULL,
  growMode = NULL,
  changeToConstant = FALSE
)
