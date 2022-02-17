### mr.ash fit/pred/const functions ###
#' @import mr.ash.alpha

weighted.mr.ash = function(X, Y, sigma2, init = NULL) {
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

  # have to scale sa2, because it is actually scaled in terms of response, I believe
  beta.init = try({rowSums(init$mr.ash.post$phi * init$mr.ash.post$m)}, silent = T)
  if (class(beta.init) == "try-error") {
    beta.init = NULL
  }
  pi.init = try({init$prior_pi})
  if (class(pi.init) == "try-error") {
    pi.init = NULL
  }
  mr.ash.fit = mr.ash(X = X_tilde, y = Y_tilde, sa2 = 1 * (2^(0.05*(0:29)) - 1)^2, max.iter = 100, beta.init = beta.init, pi = pi.init, update.sigma2 = FALSE, sigma2 = 1, standardize = FALSE, intercept = FALSE)
  mr.ash.post = get.full.posterior(mr.ash.fit)

  beta_post_1 = rowSums(mr.ash.post$phi * mr.ash.post$m)
  beta_post_2 = tcrossprod(beta_post_1)
  diag(beta_post_2) = rowSums(mr.ash.post$phi * (mr.ash.post$s2 + mr.ash.post$m^2))

  Xb_post = X %*% beta_post_1
  X_avg_b_post = sum(X_avg * beta_post_1)

  intercept = as.numeric(Y_avg - X_avg_b_post)

  # mu1 = E[int + Xb] = E[Y_avg - X_avg'b + Xb]
  mu1 = intercept + Xb_post
  # mu2 = E[(int + Xb)^2] = E[(Y_avg - X_avg'b + Xb)^2]
  #mu2 = Y_avg^2 + 2*Y_avg*(Xb_post - X_avg_b_post) + compute_X2b(X, beta_post_2, X_avg)
  mu2 = Y_avg^2 + 2*Y_avg*(Xb_post - X_avg_b_post) + apply(X_cent, MARGIN = 1, function(x) emulator::quad.form(beta_post_2, x))

  neg.elbo.mr.ash = tail(mr.ash.fit$varobj, 1)
  KL_div = neg.elbo.mr.ash -
    (.5 * length(mr.ash.fit$data$y) * log(2*pi*mr.ash.fit$sigma2)) -
    ((1 / (2 * mr.ash.fit$sigma2)) * (sum(mr.ash.fit$data$y^2) - 2*sum(mr.ash.fit$data$y * (mr.ash.fit$data$X %*% beta_post_1)) + sum(apply(mr.ash.fit$data$X, MARGIN = 1, function(x) emulator::quad.form(beta_post_2, x)))))

  return(list(mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), KL_div = KL_div, prior_pi = mr.ash.fit$pi, prior_var = mr.ash.fit$data$sa2, mr.ash.post = mr.ash.post, intercept = intercept, X_avg = X_avg, Y_avg = Y_avg))

}

fitFn.mr.ash = function(X, Y, sigma2, init) {
  return(weighted.mr.ash(X, Y, sigma2, init))
}


#' @importFrom emulator quad.form
predFn.mr.ash = function(X_new, currentFit, moment = c(1, 2)) {
  beta_post_1 = rowSums(currentFit$mr.ash.post$phi * currentFit$mr.ash.post$m)
  if (moment == 1) {
    return(currentFit$intercept + (X_new %*% beta_post_1))
  } else if (moment == 2) {
    beta_post_2 = tcrossprod(beta_post_1)
    diag(beta_post_2) = rowSums(currentFit$mr.ash.post$phi * (currentFit$mr.ash.post$s2 + currentFit$mr.ash.post$m^2))
    X_new_cent = sweep(X_new[[1]], MARGIN = 2, STATS = currentFit$X_avg, FUN = '-')
    return(currentFit$Y_avg^2 + 2*currentFit$Y_avg*(X_new %*% beta_post_1) - sum(currentFit$X_avg * beta_post_1) + apply(X_new_cent, MARGIN = 1, function(x) quad.form(beta_post_2, x)))
  } else {
    stop("`moment` must be either 1 or 2")
  }
}


constCheckFn.mr.ash = function(currentFit) {
  V = sum(currentFit$prior_pi * currentFit$prior_var) # overall variance of prior
  return(V < 1e-3)
}

# mr.ash learner
#' @export
mrAshLearner = list(
  fitFunction = fitFn.mr.ash,
  predFunction = predFn.mr.ash,
  constCheckFunction = constCheckFn.mr.ash,
  currentFit = NULL,
  X = NULL,
  X_test = NULL,
  growMode = NULL,
  changeToConstant = FALSE
)

