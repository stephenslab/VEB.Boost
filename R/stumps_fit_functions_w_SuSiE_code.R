### SuSiE stumps-related functions

calc_KL = function(mu, alpha, sigma2, V = 1) {
  p = length(mu)
  prior_weights = rep(1 / p, p)
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

optimize_V = function(tau_no_V, nu, sigma2, prior_weights, V = 1) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(nu))
  }
  lV = optim(par = log(V), fn = neg.loglik.logscale, tau_no_V = tau_no_V, nu = nu, sigma2 = sigma2, prior_weights = prior_weights, method='Brent', lower = -10, upper = 15)$par
  V = exp(lV)
  return(V)
}

#' set up a general trend filtering matrix
#' @param t vector of length n specifying locations of data points on x axis
#' @param br vector of length (p-1) specifying break points on x axis (ie where changepoints can occur)
#' By default br=t which allows breaks to occur between every data point. Note that internally duplicate elements of br are removed.
#' @param order non-negative integer indicating order of trend filtering basis (0 is changepoint basis and is the only case we test and use)
#' @keywords internal
make_tfg_matrix = function(t,br=t,order=0){
  br = unique(sort(br))
  n = length(t) # number of data points
  p = length(br) + 1 # number of bins specified by breaks
  X <- numeric(0)
  attr(X, "nrow") <- n
  attr(X, "ncol") <- p
  attr(X, "matrix.type") = "tfg_matrix"
  attr(X, "order") = order
  attr(X,"t") <- t
  attr(X,"br") <- br
  attr(X,"order_t") <- order(t)
  attr(X,"t_to_bin") <- .bincode(t,breaks = c(-Inf,br,Inf))
  attr(X,"bin_to_t") <- cumsum(hist(t, breaks = c(-Inf,br,Inf), plot=FALSE)$counts)
  attr(X,"scaled:center") <- 0
  attr(X,"scaled:scale") <- 1
  return(X)
}

is.tfg_matrix=function(X){
  ifelse(is.null(attr(X, "matrix.type")),FALSE,attr(X,"matrix.type")=="tfg_matrix")
}

# over-write make_stumps_matrix to not rely on susieR::
# also change Xtrain to be a list (allows for different lengths of breaks)
# include_linear now supports a logical vector input with the same length as ncol(X)
# for those entries that are TRUE, it includes those variables as linear terms
make_stumps_matrix = function(X, include_linear, include_stumps, Xtrain=NULL){
  if(is.null(Xtrain)){Xtrain = lapply(1:ncol(X), function(i) X[, i])}

  if (length(include_linear) == 1) { # change include_linear to be a logical vector
    include_linear = rep(include_linear, ncol(X))
  }

  xl=list() # initialize
  if(any(include_linear)){ #include X as a regular matrix first
    X_linear = X[, include_linear]
    attr(X_linear,"nrow") <- nrow(X_linear)
    attr(X_linear,"ncol") <- ncol(X_linear)
    attr(X_linear,"scaled:center") <- rep(0,ncol(X_linear))
    attr(X_linear,"scaled:scale") <- rep(1,ncol(X_linear))
    X_linear2 = X_linear^2
    attr(X_linear, "X2") <- X_linear2
    xl=c(xl,list(X_linear))
  }
  
  if (include_stumps) {
    for(i in 1:ncol(X)){
      xl= c(xl,list(make_tfg_matrix(X[,i],Xtrain[[i]])))
      }
  }

  return(xl)
}

#' @importFrom spatstat.utils revcumsum
# over-write stumps multiplication (WHY IS THERE A -1*.... ?!?!?)
#' @title Compute unscaled X \%*\% b using the special structure of trend filtering
#' @param X a tfg_matrix created by make_tfg_matrix
#' @param b a p vector of the changes at each change point
#' @return an n vector of the means at each data point
#' @keywords internal
compute_tfg_Xb = function(X,b){
  order = get_order(X)
  for(i in 1:(order+1)){
    #b = rev(cumsum(rev(b))) # computes mean in each bin
    b = spatstat.utils::revcumsum(b) # faster than rev(cumsum(rev(b)))
  }
  return(b[attr(X,"t_to_bin")]) #  maps bin means to a mean for each datapoint
}

#' @title Compute t(X) \%*\% y using the special structure of trend filtering
#' @param X a tfg_matrix created by make_tfg_matrix
#' @param y an n vector of data
#' @return a p vector
#' @keywords internal
compute_tfg_Xty = function(X,y){
  order = get_order(X)
  y = y[attr(X,"order_t")] # sort y according to increasing t
  for (i in 1:(order+1)){
    y = cumsum(y)
  }
  return(y[attr(X,"bin_to_t")])
}

# computes (X - X_avg)^2 %*% b
compute_X2b = function(X, b, X_avg = 0) {
  if (is.list(X)) {
    n_var = sapply(X, get_ncol) # number of variables for each element of X
    b_split = split_vector(b, n_var) # split b into a list of vectors
    X_avg_split = split_vector(X_avg, n_var)
    X2b = mapply(compute_X2b, X, b_split, X_avg_split, SIMPLIFY = FALSE) # apply compute_X2b to elements of lists X, b_split
    return(Reduce(`+`, X2b)) # add the results up
  } else {
    if (is.tfg_matrix(X)) {
      # X is boolean matrix, so X^2 = X
      return(compute_Xb(X, b) - 2*compute_Xb(X, b*X_avg) + sum(X_avg^2 * b))
    } else {
      return(compute_Xb(attr(X, "X2"), b) - 2*compute_Xb(X, b*X_avg) + sum(X_avg^2 * b))
    }
  }
}

# computes t((X - X_avg)^2) %*% y
compute_X2ty = function(X, y, X_avg = 0) {
  if (is.list(X)) {
    X_avg_split = split_vector(X_avg, sapply(X, get_ncol))
    y_list = rep(list(y), length(X_avg_split))
    return(unlist(mapply(compute_X2ty, X, y_list, X_avg_split)))
  } else {
    if (is.tfg_matrix(X)) {
      # X is boolean matrix, so X^2 = X
      return(as.numeric(compute_Xty(X, y)) * (1 - 2*X_avg) + (X_avg^2 * sum(y)))
    } else {
      return(as.numeric(compute_Xty(attr(X, "X2"), y) - 2*compute_Xty(X, y)*X_avg + (X_avg^2 * sum(y))))
    }
  }
}


# weighted SER function, linear terms + stumps
# X is a list, first element corresponds to linear, others are stumps for variables
weighted_SER = function(X, Y, sigma2, init = list(V = NULL)) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(Y))
  }

  inv_sigma2 = 1 / sigma2
  sum_inv_sigma2 = sum(inv_sigma2)
  w = inv_sigma2 / sum_inv_sigma2
  p = get_ncol(X)
  prior_weights = rep(1 / p, p)
  Y_avg = sum(Y * w)
  Y_cent = Y - Y_avg
  X_avg = compute_Xty(X, w) # vector of weighted avg of columns of X

  tau_no_V = compute_X2ty(X, inv_sigma2, X_avg)
  nu = compute_Xty(X, Y_cent / sigma2) - (X_avg * sum(Y_cent / sigma2))

  # optim method, seems to be slower than EM method
  V = ifelse(is.null(init$V), 1, init$V)
  V = optimize_V(tau_no_V, nu, sigma2, prior_weights, V)

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

  KL_div = calc_KL(mu, alpha, sigma2_post, V)

  return(list(mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), KL_div = KL_div, alpha = alpha, mu = mu, sigma2_post = sigma2_post, intercept = intercept, V = V, X_avg = X_avg, Y_avg = Y_avg))
}

fitFnSusieStumps = function(X, Y, sigma2, init) {
  return(weighted_SER(X, Y, sigma2, init))
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
