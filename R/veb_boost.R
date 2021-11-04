#' Performs VEB-Boosting
#'
#' Solves the VEB-Boost regression problem using the supplied inputs
#'
#' @details
#'
#' Given a pre-specified arithmetic tree structure \deqn{T(\mu_1, \dots, \mu_L)}, \deqn{\mu_l := h_l(\beta_l)}, 
#' priors \deqn{\beta_l \sim g_l(\cdot)}, and inputs for the response, VEB-Boosting is performed.
#'
#' A cyclic CAVI scheme is used, where we cycle over the leaf nodes and update the approxiomation
#' to the posterior distribution at each node in turn.
#'
#' We start with the arithmetic tree structure \deqn{T(\mu_1, \dots, \mu_L) = \sum_{i=1}^k \prod_{j=1}^{d_k} \mu_{i, j}}
#'
#'
#' @param X is a predictor object to be used.
#' This object can take any form, so long as the user-supplied \code{fitFunctions} and \code{predFunctions} know how to use them.
#'
#' @param Y is a numeric response. For all but the 'aft.loglogistic' family, this should be an n-vector.
#' For the 'aft.loglogistic' family, this should be an n x 2 matrix, with the first column being the left end-point of survival time
#' and the second column being the right end-point of survival time (used for interval-censored observations)..
#' In the case of left-censored data, the left end-point should be 'NA', and in the case of right-censored data, the right end-point should be 'NA'.
#' If the observation is uncensored, both end-points should be equal to the observed survival time.
#' 
#' @param X_test is an optional predictor object to be used as the testing data. Posterior mean response is saved in the output's field \code{$pred_mu1}.
#' Alternatively, after running \code{veb_boost}, the user can call the \code{$predict} method on the output, with \code{X_test} as the first
#' argument, and \code{1} as the second argument.
#'
#' @param fitFunctions is either a single fitting function, or a list of length \code{k} of fitting functions to be used in
#' each term on the sum of nodes
#'
#' @param predFunctions is either a single prediction function, or a list of length \code{k} of prediction functions to be used in
#' each term of the sum of nodes
#'
#' @param constCheckFunctions is either a single constant check function, or a list of length \code{k} of constant check functions
#' to be used in each term of the sum of nodes
#' 
#' @param extrapolateFunctions is either a single fit extrapolation function, or a list of length \code{k} of fit extrapolation functions
#' to be used in each term of the sum of nodes
#'
#' @param growTree is a logical for if we should grow the tree after convergence (TRUE), or only use the initial tree
#' structure (FALSE)
#'
#' @param k is an integer for how many terms are in the sum of nodes
#'
#' @param d is either an integer, or an integer vector of length \code{k} for the multiplicative depth of each of the k terms
#' NOTE: This can be dangerous. For example, if the fit starts out too large, then entire branhces will be fit to be exactly
#' zero. When this happens, we end up dividing by 0 in places, and this results in NAs, -Inf, etc. USE AT YOUR OWN RISK
#'
#' @param growMode specifies how we grow the tree, either splitting nodes with addition, multiplication, or both
#' If \code{+*}, we grow mu_0 -> (mu_0 * mu_2) + mu_1
#' If \code{+}, we grow mu_0 -> (mu_0 + mu_1)
#' If \code{*}, we grow mu_0 -> (mu_0 * mu_1) (NOTE: Not recommended if we start with \code{k = 1})
#' 
#' @param sigma2 is a scalar/n-vector specifying a fixed residual variance. If not NULL, then the residual variance will be
#' fixed to this value/vector. If NULL, then it will be initialized and updated automatically.
#' This should be left as NULL unless you really know what you're doing. For safety, this can only be not NULL if family is gaussian
#' 
#' @param addMrAsh a logical flag for if a Mr.Ash fit should be added to the full tree to "regress out" linear effects
#' 
#' @param R is a correlation matrix for any random effect to be added to the VEB-Boost tree (if NULL, no random intercept is included)
#'
#' @param changeToConstant is a flag for if, when the fit is found to be basically constant, if we should actually change
#' the fitting function of that node to fit exactly a constant value
#'
#' @param family is what family the response is
#'
#' @param tol is a positive scalar specifying the level of convergence to be used
#' 
#' @param exposure is a scalar or a vector used for the Poisson regression, NB, and AFT cases. For Poisson, we assume that the response
#' \deqn{Y_i \sim Pois(c_i \lambda_i)} for given exposure \deqn{c_i}, and we model \deqn{\lambda_i}
#' For AFT, exposure is 1 for non-censored observations, and 0 for right-censored observations
#'
#' @param verbose is a logical flag specifying whether we should report convergence information as we go
#' 
#' @param maxit is the maximum number of iterations for each version of the VEB-Boost tree
#' 
#' @param backfit is a logical. If TRUE, then after the algorithm is done, it'll run through once more with the
#' current tree and fit to convergence. Useful when, e.g. maxit = 1. 
#'
#' @param mc.cores is the number of cores to use in mclapply, only used in family == "multinomial", and only
#' supported on UNIX systems, where mclapply works. NOT CURRENTLY SUPPORTED
#'
#'
#' @return A \code{VEB_Boost_Node} object with the fit
#' 
#' 
#' @examples
#' 
#' set.seed(1)
#' n = 1000
#' p = 1000
#' X = matrix(runif(n * p), nrow = n, ncol = p)
#' Y = rnorm(n, 5*sin(3*X[, 1]) + 2*(X[, 2]^2) + 3*X[, 3]*X[, 4])
#' 
#' For input X and list \code{fit} returned from fitFn (encoding the variational posterior for b), computes either the 1st or 2nd 
#' predFn = function(X, fit, moment = c(1, 2)) {
#'   # posterior moments of the response (depending on if \code{moment} is 1 or 2)
#'   if (moment == 1) {
#'     res = E_fit[Xb]
#'   } else {
#'     res = E_fit[(Xb)^2]
#'   }
#'   return(res)
#' }
#' 
#' For a given prior g(b), a function that approximates the posterior of b, q(b), using Variational Inference
#' fitFn = function(X, Y, sigma2, init) {
#'   fit = list(whatever is needed to encode the variational posterior)
#'   KL_div = D_KL(q || g) # KL divergence from variational posterior to prior
#'   mu1 = predFn(X, fit, 1)
#'   mu2 = predFn(X, fit, 2)
#'   # add mu1, mu2, and KL_div to the fit, MUST BE CALLED $mu1, $mu1, and $KL_div
#'   fit$mu1 = mu1
#'   fit$mu2 = mu2
#'   fit$KL_div = KL_div
#'   return(fit)
#' }
#' 
#' For a given \code{fit}, returns TRUE if the variational posterior is close enough to a constant, else FALSE
#' constCheckFn = function(fit) {
#'   if (fit is close to constant) {
#'     return(TRUE)
#'   } else {
#'     return(FALSE)
#'   }
#' }
#' veb.fit = veb_boost(X, Y, fitFn, predFn, constCheckFn, family = "gaussian")
#'
#' @export
#'

veb_boost = function(X, Y, X_test = NULL, fitFunctions, predFunctions, constCheckFunctions,
                     growTree = TRUE, k = 1, d = 1, growMode = c("+*", "+", "*"), sigma2 = NULL, addMrAsh = FALSE, R = NULL,
                     changeToConstant = FALSE, family = c("gaussian", "binomial", "multinomial", "negative.binomial", "poisson.log1pexp", "poisson.exp", "aft.loglogistic", "ordinal.logistic"),
                     exposure = NULL, tol = nrow(as.matrix(Y)) / 10000, verbose = TRUE, maxit = Inf, backfit = FALSE, mc.cores = 1) {

  ### Check Inputs ###
  if (!(is.vector(Y) || ((family == 'aft.loglogistic') && (is_valid_matrix(Y)) && (ncol(Y) == 2)))) {
    stop("'Y' must be a vector (or a matrix with 2 columns for the 'aft.loglogistic' case)")
  }
  if (!is.null(sigma2)) {
    if (family != "gaussian") {
      stop("'sigma2' can only be not NULL if 'family' is set to 'gaussian'")
    }
    if (!(length(sigma2) %in% c(1, nrow(as.matrix(Y))))) {
      stop("'sigma2' must have length 1 or length of Y")
    }
    if (any(sigma2 <= 0)) {
      stop("'sigma2' must contain only positive values")
    }
  }
  # Logical Flags
  if (!(growTree %in% c(TRUE, FALSE))) {
    stop("'growTree' must be either TRUE or FALSE")
  }
  if (!(verbose %in% c(TRUE, FALSE))) {
    stop("'verbose' must be either TRUE or FALSE")
  }
  if (!(changeToConstant %in% c(TRUE, FALSE))) {
    stop("'changeToConstant' must be either TRUE or FALSE")
  }
  if (!(backfit %in% c(TRUE, FALSE))) {
    stop("'backfit' must be either TRUE or FALSE")
  }
  # Scalars
  if ((k < 1) || (k %% 1 != 0)) {
    stop("'k' must be a positive whole number")
  }
  if ((d < 1) || (d %% 1 != 0)) {
    stop("'d' must be a positive whole number")
  }
  # Functions
  if (class(fitFunctions) == "function") {
    fitFunctions = list(fitFunctions)
  }
  if ((class(fitFunctions) != "list") || !(length(fitFunctions) %in% c(1, k)) || (class(fitFunctions[[1]]) != "function")) {
    stop("'fitFunctions' must either be a function, or a list of functions of length 1 or k")
  }
  if (class(predFunctions) == "function") {
    predFunctions = list(predFunctions)
  }
  if ((class(predFunctions) != "list") || !(length(predFunctions) %in% c(1, k)) || (class(predFunctions[[1]]) != "function")) {
    stop("'predFunctions' must either be a function, or a list of functions of length 1 or k")
  }
  if (class(constCheckFunctions) == "function") {
    constCheckFunctions = list(constCheckFunctions)
  }
  if ((class(constCheckFunctions) != "list") || !(length(constCheckFunctions) %in% c(1, k)) || (class(constCheckFunctions[[1]]) != "function")) {
    stop("'constCheckFunctions' must either be a function, or a list of functions of length 1 or k")
  }
  # Grow Mode
  growMode = match.arg(growMode)
  # Family
  family = match.arg(family)
  # Exposure
  if (!(grepl("poisson", family, ignore.case = TRUE) | family == "negative.binomial")) {
    if (!is.null(exposure)) {
      warning("Argument 'exposure' is only used in the Poisson, or negative binomial cases, ignoring supplied input")
    }
  } else {
    if (is.null(exposure)) {
      warning("Argument 'exposure' must be supplied in the Poisson, or negative binomial cases, setting to 1")
      exposure = rep(1, length(Y))
    }
  }
  if (!is.null(exposure)) {
    if (any(exposure < 0)) {
      stop("When using the poisson on NB families, 'exposure' must be positive")
    }
  }
  # RE correlation matrix
  if (!is.null(R)) {
    # if (!is_valid_matrix(R) || !Matrix::isSymmetric(R) || (max(abs(R)) > 1) || any(diag(R) != 1)) {
    if (!is_valid_matrix(R) || !Matrix::isSymmetric(R) || (max(abs(R)) > 1)) {
      stop("'R' must be a correlation matrix")
    }
    if (is.null(attr(R, 'svd'))) {
      stop("'R' must have an attribute names 'svd' that contains the (reduced) SVD of R (see, e.g. sparsesvs::sparsesvd)")
    }
  }
  # Tolerance
  if (!is.numeric(tol) || (length(tol) > 1) || (tol <= 0)) {
    stop("'tol' must be a positive number")
  }
  # Response
  if ((family == "gaussian") && !is.numeric(Y)) {
    stop("'Y' must be numeric when using gaussian response")
  }
  if ((family == "binomial") && any(!(Y %in% c(0, 1)))) {
    stop("'Y' must be 0/1 when using binomial response")
  }
  if ((grepl("poisson", family, ignore.case = TRUE) | family == "negative.binomial") && (any(Y < 0) || any(Y %% 1 != 0))) {
    stop("'Y' must be positive integers when using Poisson or negative binomial response")
  }
  if ((family == "aft.loglogistic")) {
    if (any(Y <= 0, na.rm = TRUE)) {
      stop("Non-NA 'Y' must be positive when using the AFT model")
    }
    if (any(Y[, 2] - Y[, 1] < 0, na.rm = TRUE)) {
      stop("The first column of 'Y' must be <= the second column of 'Y'")
    }
    if (any(rowSums(is.na(Y)) == 2)) {
      stop("'Y' cannot have rows that are all NA")
    }
  }
  if (family == "ordinal.logistic") {
    if (any(!(Y %in% (1:max(Y))))) {
      stop("'Y' must contain only integer values from 1 to max(Y)")
    }
  }
  # mc.cores
  if (mc.cores != 1) {
    message("Parallel multinomial updates not currently supported, setting 'mc.cores' to 1")
    mc.cores = 1
  }
  # maxit
  if (maxit < 1) {
    stop("'maxit' must be a positive integer")
  }
  maxit = ceiling(maxit)
  

  # get censoring types for AFT model
  if (family == "aft.loglogistic") {
    attr(Y, 'left_cens') = (is.na(Y[, 1]) & !is.na(Y[, 2]))
    attr(Y, 'right_cens') = (!is.na(Y[, 1]) & is.na(Y[, 2]))
    attr(Y, 'int_cens') = (!is.na(rowSums(Y)) & (Y[, 1] < Y[, 2]))
    attr(Y, 'not_cens') = (!is.na(rowSums(Y)) & (Y[, 1] == Y[, 2]))
    attr(Y, 'log_Y') = log(Y)
  }
  if (family == "ordinal.logistic") {
    attr(Y, 'left_cens') = (Y == 1)
    attr(Y, 'right_cens') = (Y == max(Y))
    attr(Y, 'int_cens') = !(attr(Y, 'left_cens') | attr(Y, 'right_cens'))
  }
  
  ### Run Gaussian/Binomial Cases ###
  if (family %in% c("gaussian", "binomial", "negative.binomial", "poisson.log1pexp", "poisson.exp", "aft.loglogistic", "ordinal.logistic")) {
    # initialize tree
    learner = initialize_veb_boost_tree(X = X, Y = Y, k = k, d = d, fitFunctions = fitFunctions, predFunctions = predFunctions, constCheckFunctions = constCheckFunctions, 
                                        addMrAsh = addMrAsh, R = R, family = family, exposure = exposure)
    if (family == "gaussian") {
      learner$sigma2 = if (is.null(sigma2)) var(Y) else sigma2
    } else if (family == "aft.loglogistic") {
      learner$sigma2 = var(rowMeans(attr(Y, 'log_Y'), na.rm = TRUE))
    } else if (family == "ordinal.logistic") {
      learner$cutpoints = seq(from = -max(Y), to = max(Y), length.out = max(Y) - 1)
      attr(learner$cutpoints, "log_Y") = matrix(NA, nrow = length(Y), ncol = 2)
      attr(learner$cutpoints, "log_Y")[attr(Y, 'left_cens'), 2] = learner$cutpoints[1]
      for (k in 2:(length(learner$cutpoints))) {
        attr(learner$cutpoints, "log_Y")[which(Y == k), 1] = learner$cutpoints[k-1]
        attr(learner$cutpoints, "log_Y")[which(Y == k), 2] = learner$cutpoints[k]
      }
      attr(learner$cutpoints, "log_Y")[attr(Y, 'right_cens'), 1] = tail(learner$cutpoints, 1)
    }
    update_sigma2 = ifelse(is.null(sigma2), family %in% c("gaussian", "aft.loglogistic", "ordinal.logistic"), FALSE) # only update variance if Gaussian , or AFT model

    # converge fit once
    learner$convergeFit(tol, update_sigma2 = update_sigma2, verbose = FALSE, maxit = maxit)

    # if growing tree, continue w/ growing tree & fitting to convergence
    if (growTree) {
      if (verbose) {
        cat(paste("ELBO: ", round(learner$ELBO, 3), sep = ""))
        cat("\n")
      }
      
      learner$convergeFitAll(tol = tol, update_sigma2 = update_sigma2, growMode = growMode, changeToConstant = changeToConstant, verbose = FALSE, maxit = maxit)

      # while ((tail(tail(learner$ELBO_progress, 1)[[1]], 1) - tail(tail(learner$ELBO_progress, 2)[[1]], 1) > tol) && 
      #        (length(Traverse(learner, filterFun = function(node) node$isLeaf & !node$isLocked)) > 0) && (tail(tail(learner$ELBO_progress, 1)[[1]], 1) != Inf)) {
      while ((mean(diff(sapply(tail(learner$ELBO_progress, 3), tail, 1))) > 10*tol) && 
             (length(Traverse(learner, filterFun = function(node) node$isLeaf & !node$isLocked)) > 0) && (tail(tail(learner$ELBO_progress, 1)[[1]], 1) != Inf)) {
        if (verbose) {
          cat(paste("ELBO: ", round(learner$ELBO, 3), sep = ""))
          cat("\n")
        }
        learner$convergeFitAll(tol = tol, update_sigma2 = update_sigma2, growMode = growMode, changeToConstant = changeToConstant, verbose = FALSE, maxit = maxit)
      }
      
    }
    
    learner$lockLearners(growMode = growMode, changeToConstant = changeToConstant) # change to constant at the end if needed
    if (backfit) { # converge again, if backfitting
      learner$convergeFit(tol, update_sigma2 = update_sigma2, verbose = FALSE, maxit = Inf)
    }
    
    if (!is.null(X_test)) {
      learner$predict(X_test, 1)
    }
    return(learner)
    
  } else { # else, multinomial case
    classes = sort(unique(Y))
    learner_multiclass = VEBBoostMultiClassLearner$new()
    learner_multiclass$Y = Y
    learner_multiclass$mc.cores = mc.cores
    learner_multiclass$classes = classes
    learner_multiclass$X = X

    learnerList = list()
    for (j in 1:length(classes)) { # add learners to learnerList
      learner = initialize_veb_boost_tree(X = NULL, Y = 1 * (Y == classes[j]), k = k, d = d, fitFunctions = fitFunctions, predFunctions = predFunctions,
                                constCheckFunctions = constCheckFunctions, family = "binomial")
      learnerList = c(learnerList, learner)
    }
    names(learnerList) = paste0("learner", 0:(length(classes) - 1), sep = "")
    learner_multiclass$addLearners(learnerList) # add learner list to multi-class clearner

    # converge fit once
    learner_multiclass$convergeFit(tol = tol)

    # if growing tree, continue w/ growing tree & fitting to convergence
    if (growTree) {
      if (verbose) {
        cat(paste("ELBO: ", round(learner_multiclass$ELBO, 3), sep = ""))
        cat("\n")
      }
      
      learner_multiclass = learner_multiclass$convergeFitAll(tol = tol, growMode = growMode, changeToConstant = changeToConstant)

      while ((abs(tail(tail(learner_multiclass$ELBO_progress, 1)[[1]], 1) - tail(tail(learner_multiclass$ELBO_progress, 2)[[1]], 1)) > tol) &&
             (sum(sapply(learner_multiclass$learners, function(x) length(Traverse(x, filterFun = function(node) node$isLeaf & !node$isLocked)))) > 0)) {
        if (verbose) {
          cat(paste("ELBO: ", round(learner_multiclass$ELBO, 3), sep = ""))
          cat("\n")
        }
        learner_multiclass$convergeFitAll(tol = tol, growMode = growMode, changeToConstant = changeToConstant)
      }
    }
    
    if (!is.null(X_test)) {
      learner_multiclass$predict(X_test, 1)
    }

    return(learner_multiclass)
  }

}


#' Wrapper for using VEB-Boost with the SER prior w/ stumps and/or linear terms
#' 
#' @details
#'
#' This function performs VEB-Boosting, where the prior to be used is the SER prior, and our predictors are either
#' i) the linear terms of X;
#' ii) the stumps made from the columns of X; or
#' iii) Both (i) and (ii)
#' 
#' @param X An (n x p) numeric matrix to be used as the predictors (currently, this wrapper forces all nodes to use the same X)
#' 
#' @param Y is a numeric vector response
#' 
#' @param X_test is an optional (m X p) matrix to be used as the testing data. Posterior mean response is saved in the output's field \code{$pred_mu1}
#' 
#' @param include_linear is a logical of length 1 or p specifying which columns of X we should include as linear terms.
#' If the length is 1, this value gets recycled for all columns of X.
#' If NULL is supplied, then all valid linear terms are used.
#' 
#' @param include_stumps is a logical of length 1 or p specifying which columns of X we should include as stump terms
#' If the length is 1, this value gets recycled for all columns of X.
#' If NULL is supplied, then all valid stumps terms are used.
#' 
#' @param num_cuts is a whole number of length 1 or p specifying how many cuts to make when making the stumps terms.
#' If the length is 1, this value gets recycled for all columns of X.
#' For entries corresponding to the indices where \code{include_stumps} is FALSE, these values are ignored.
#' We use the quantiles from each predictor when making the stumps splits, using \code{num_cuts} of them.
#' If \code{num_cuts = Inf}, then all values of the variables are used as split points.
#' 
#' @param use_quants is a logical for if the cut-points should be based off of the quantiles (`use_quants = TRUE`), or if the
#' cut points should be evenly spaced in the range of the variable (`use_quants = FALSE`).
#' 
#' @param scale_X is a string for if/how the columns of X should be scaled.
#' 'sd' scales by the standard deviations of the variables.
#' 'max' scales by the maximum absolute value (so variables are on the [-1, +1] scale).
#' 'NA' performs no scaling.
#' 
#' @param max_log_prior_var is a scalar for the maximum that the estimated log-prior variance for each weak learner can be.
#' The idea is that setting this to be small limits the "size" of each weak learner, similar-ish to the learning rate in boosting.
#' The maximum allowed value is 35, which essentially allows the unrestricted MLE to be estimated.
#' The minimum allowed value is -5.
#' 
#' @param lin_prior_prob is a number between 0 and 1 that gives the prior probability that the effect variable is a linear term.
#' This means that (1 - lin_prior_prob) is the probability that the effect variable is a stump term.
#' Within linear terms and stump terms, all variables have the same prior variance.
#' 
#' @param ... Other arguments to be passed to \code{\link{veb_boost}}
#' 
#' @return A \code{VEB_Boost_Node} object with the fit
#'
#' @examples 
#' set.seed(1)
#' n = 1000
#' p = 1000
#' X = matrix(runif(n * p), nrow = n, ncol = p)
#' Y = rnorm(n, 5*sin(3*X[, 1]) + 2*(X[, 2]^2) + 3*X[, 3]*X[, 4])
#' veb.stumps.fit = veb_boost_stumps(X, Y, include_linear = TRUE, family = "gaussian")
#' 
#' @importFrom sparseMatrixStats colRanges
#' @importFrom sparseMatrixStats colQuantiles
#' 
#' @useDynLib VEB.Boost, .registration=TRUE
#' 
#' @export
#' 

veb_boost_stumps = function(X, Y, X_test = NULL, include_linear = NULL, include_stumps = NULL, 
                            num_cuts = ceiling(min(length(Y) / 5, max(100, sqrt(length(Y))))), 
                            use_quants = TRUE, scale_X = c("sd", "max", "NA"), 
                            max_log_prior_var = 0, lin_prior_prob = 0.5, nthreads = ceiling(parallel::detectCores(logical = FALSE) / 2), ...) {
  # # Set higher priority for this process
  # psnice_init = tools::psnice()
  # on.exit({tools::psnice(tools::psnice(value = psnice_init))})
  # tools::psnice(value = 0 - 10*(Sys.info()["sysname"] == "Windows"))

  ### Check Inputs ###
  # Check X
  if (!is_valid_matrix(X)) {
    stop("'X' must be a numeric matrix")
  }
  if (!is.null(X_test) && !is_valid_matrix(X_test)) {
    stop("'X_test' must be a numeric matrix")
  }
  p = ncol(X)
  # Check logicals
  if (!(is.null(include_linear) || all(include_linear %in% c(TRUE, FALSE)))) {
    stop("'include_linear' must be either TRUE or FALSE, or be NULL")
  }
  if (!(is.null(include_linear) || (length(include_linear) %in% c(1, p)))) {
    stop("'include_linear' must have length 1 or ncol(X), or be NULL")
  }
  if (!(is.null(include_stumps) || all(include_stumps %in% c(TRUE, FALSE)))) {
    stop("'include_stumps' must be either TRUE or FALSE, or be NULL")
  }
  if (!(is.null(include_stumps) || (length(include_stumps) %in% c(1, p)))) {
    stop("'include_stumps' must have length 1 or ncol(X), or be NULL")
  }
  # Make sure at least one include is true
  if (!is.null(include_linear) && !is.null(include_stumps) && all(include_linear == FALSE) && all(include_stumps == FALSE)) {
    stop("At least one of 'include_linear' or 'include_stumps' must be TRUE or NULL")
  }
  # Make sure num_cuts is a positive whole number
  if (any(num_cuts < 1) || any(num_cuts[is.finite(num_cuts)] %% 1 != 0)) {
    stop("'num_cuts' must be a positive whole number or Inf")
  }
  # if (!(length(num_cuts) %in% c(1, p))) {
  #   stop("'num_cuts' must have length 1 or ncol(X)")
  # }
  if (length(num_cuts) != 1) {
    stop("'num_cuts' must have length 1")
  }
  if (!(use_quants %in% c(TRUE, FALSE))) {
    stop("'use_quants' must be either TRUE or FALSE")
  }
  scale_X = match.arg(scale_X)
  # Check max_log_prior_var
  if (!is.numeric(max_log_prior_var) || (length(max_log_prior_var) != 1)) {
    stop("'max_log_prior_var' must be a single number")
  }
  if ((max_log_prior_var < -5) | (max_log_prior_var > 35)) {
    warning("Restricting 'max_log_prior_var' to be in the range [-5, 35]")
    max_log_prior_var = min(35, max(-5, max_log_prior_var))
  }
  # Check lin_prior_prob
  if (!is.numeric(lin_prior_prob) || (length(lin_prior_prob) != 1)) {
    stop("'lin_prior_prob' must be a single number")
  }
  if ((lin_prior_prob < 0) | (lin_prior_prob > 1)) {
    warning("Restricting 'lin_prior_prob' to be in the range [0, 1]")
    lin_prior_prob = min(1, max(0, lin_prior_prob))
  }
  
  # scale X if needed
  if (scale_X %in% c("max", "sd")) {
    if (scale_X == "max") {
      X_scale_factors = sparseMatrixStats::colMaxs(abs(X))
    } else {
      X_scale_factors = sparseMatrixStats::colSds(X)
    }
    if (inherits(X, "CsparseMatrix")) {
      X@x = X@x / rep.int(X_scale_factors, diff(X@p))
    } else {
      X = sweep(X, 2, X_scale_factors, '/')
    }
    if (!is.null(X_test)) {
      if (inherits(X_test, "CsparseMatrix")) {
        X_test@x = X_test@x / rep.int(X_scale_factors, diff(X_test@p))
      } else {
        X_test = sweep(X_test, 2, X_scale_factors, '/')
      }
    }
  }
  # get which columns to include as linear and stumps terms, if not provided
  if (is.null(include_linear) | is.null(include_stumps)) {
    n_unique = apply(X, 2, function(x) length(unique(x)))
    if (is.null(include_linear)) {
      include_linear = (n_unique > 1)
    }
    if (is.null(include_stumps)) {
      include_stumps = (n_unique > 2)
    }
  }
  # Make stumps matrix
  if (is.infinite(num_cuts) | all(!include_stumps)) {
    cuts = NULL
  } else {
    # if (length(num_cuts) == 1) {
    #   num_cuts = rep(num_cuts, p)
    # }
    # cuts = rep(list(NULL), p)
    # for (j in 1:p) {
    #   if (is.finite(num_cuts[j])) {
    #     if (use_quants) {
    #       cuts[[j]] = quantile(X[, j], probs = qbeta(seq(from = 0, to = 1, length.out = num_cuts[j] + 2), .5, .5))[-c(1, num_cuts[j] + 2)]
    #     } else {
    #       cuts[[j]] = seq(from = min(X[, j]), to = max(X[, j]), length.out = num_cuts[j] + 2)[-c(1, num_cuts[j] + 2)]
    #     }
    #   }
    # }
    ppts = ppoints(num_cuts)
    if (use_quants) {
      col_quants = sparseMatrixStats::colQuantiles(X, probs = ppts)
      # cuts = lapply(1:p, function(j) unique(col_quants[j, ]))
      cuts = asplit(col_quants, 1)
    } else {
      col_ranges = sparseMatrixStats::colRanges(X)
      cuts = lapply(1:p, function(j) col_ranges[j, 1]*(1 - ppts) + col_ranges[j, 2]*ppts)
    }
  }
  if (inherits(X, "dgCMatrix")) {
    X_stumps = make_stumps_matrix_sp_cpp(X, 1*include_linear, 1*include_stumps, cuts, nthreads)
  } else {
    X_stumps = make_stumps_matrix_cpp(X, 1*include_linear, 1*include_stumps, cuts, nthreads)
  }  
  
  # set up testing data, if any
  X_test_stumps = NULL
  if (!is.null(X_test)) {
    # re-define cuts, for case where cuts = NULL, so that we use training data for splits, not test data
    if (is.null(cuts)) {
      brs = lapply(X_stumps, function(x) attr(x, 'br'))
      if (any(include_linear)) {
        brs = brs[-1]
      }
      cuts = vector('list', p)
      cuts[which(include_stumps)] = brs
    }
    if (inherits(X_test, "dgCMatrix")) {
      X_test_stumps = make_stumps_matrix_sp_cpp(X_test, 1*include_linear, 1*include_stumps, cuts, nthreads)
    } else {
      X_test_stumps = make_stumps_matrix_cpp(X_test, 1*include_linear, 1*include_stumps, cuts, nthreads)
    } 
  }
  
  fitFnSusieStumps_maxlV = function(X, Y, sigma2, init) {
    return(weighted_SER_cpp(X, Y, sigma2, init, max_log_prior_var, lin_prior_prob))
  }

  # Run
  veb.fit = veb_boost(X = X_stumps, Y = Y, X_test = X_test_stumps, 
                      fitFunctions = fitFnSusieStumps_maxlV, predFunctions = predFnSusieStumps_cpp, constCheckFunctions = constCheckFnSusieStumps, 
                      ...)
  
  return(veb.fit)
}


#' #' Wrapper for using VEB-Boost with the FLAM prior w/ stumps
#' #' 
#' #' @details
#' #'
#' #' This function performs VEB-Boosting, where the prior to be used is the FLAM prior, and our predictors are the stumps made from the columns of X
#' #' 
#' #' @param X An (n x p) numeric matrix to be used as the predictors (currently, this wrapper forces all nodes to use the same X)
#' #' 
#' #' @param Y is a numeric vector response
#' #' 
#' #' @param X_test is an optional (m X p) matrix to be used as the testing data. Posterior mean response is saved in the output's field \code{$pred_mu1}
#' #' 
#' #' @param num_cuts is a whole number of length 1 or p specifying how many cuts to make when making the stumps terms.
#' #' If the length is 1, this value gets recycled for all columns of X.
#' #' For entries corresponding to the indices where \code{include_stumps} is FALSE, these values are ignored.
#' #' We use the quantiles from each predictor when making the stumps splits, using \code{num_cuts} of them.
#' #' If \code{num_cuts = Inf}, then all values of the variables are used as split points.
#' #' 
#' #' @param use_quants is a logical for if the cut-points should be based off of the quantiles (`use_quants = TRUE`), or if the
#' #' cut points should be evenly spaced in the range of the variable (`use_quants = FALSE`).
#' #' 
#' #' 
#' #' @param ... Other arguments to be passed to \code{\link{veb_boost}}
#' #' 
#' #' @return A \code{VEB_Boost_Node} object with the fit
#' #'
#' #' @examples 
#' #' set.seed(1)
#' #' n = 1000
#' #' p = 1000
#' #' X = matrix(runif(n * p), nrow = n, ncol = p)
#' #' Y = rnorm(n, 5*sin(3*X[, 1]) + 2*(X[, 2]^2) + 3*X[, 3]*X[, 4])
#' #' veb.flam.fit = veb_boost_flam(X, Y, family = "gaussian")
#' #' 
#' #' @importFrom sparseMatrixStats colRanges
#' #' @importFrom sparseMatrixStats colQuantiles
#' #' 
#' #' 
#' #' 
#' 
#' veb_boost_flam = function(X, Y, X_test = NULL, 
#'                             num_cuts = ceiling(min(length(Y) / 5, max(100, sqrt(length(Y))))), 
#'                             use_quants = TRUE, ...) {
#'   ### Check Inputs ###
#'   # Check X
#'   if (!is_valid_matrix(X)) {
#'     stop("'X' must be a numeric matrix")
#'   }
#'   if (!is.null(X_test) && !is_valid_matrix(X_test)) {
#'     stop("'X_test' must be a numeric matrix")
#'   }
#'   p = ncol(X)
#'   # Make sure num_cuts is a positive whole number
#'   if (any(num_cuts < 1) || any(num_cuts[is.finite(num_cuts)] %% 1 != 0)) {
#'     stop("'num_cuts' must be a positive whole number or Inf")
#'   }
#'   if (length(num_cuts) != 1) {
#'     stop("'num_cuts' must have length 1")
#'   }
#'   if (!(use_quants %in% c(TRUE, FALSE))) {
#'     stop("'use_quants' must be either TRUE or FALSE")
#'   }
#' 
#'   # Make stumps matrix
#'   if (is.infinite(num_cuts)) {
#'     cuts = NULL
#'   } else {
#'     ppts = ppoints(num_cuts)
#'     if (use_quants) {
#'       col_quants = sparseMatrixStats::colQuantiles(X, probs = ppts)
#'       cuts = asplit(col_quants, 1)
#'     } else {
#'       col_ranges = sparseMatrixStats::colRanges(X)
#'       cuts = lapply(1:p, function(j) col_ranges[j, 1]*(1 - ppts) + col_ranges[j, 2]*ppts)
#'     }
#'   }
#'   X_stumps = make_stumps_matrix(X, FALSE, TRUE, cuts)
#'   
#'   
#'   # set up testing data, if any
#'   X_test_stumps = NULL
#'   if (!is.null(X_test)) {
#'     # re-define cuts, for case where cuts = NULL, so that we use training data for splits, not test data
#'     if (is.null(cuts)) {
#'       cuts = lapply(X_stumps, function(x) attr(x, 'br'))
#'       if (any(include_linear)) {
#'         cuts = cuts[-1]
#'       }
#'     }
#'     X_test_stumps = make_stumps_matrix(X_test, FALSE, TRUE, cuts)
#'   }
#'   
#'   # Run
#'   veb.fit = veb_boost(X = X_stumps, Y = Y, X_test = X_test_stumps, 
#'                       fitFunctions = fitFnFLAM, predFunctions = predFnFLAM, constCheckFunctions = constCheckFnFLAM,
#'                       ...)
#'   
#'   return(veb.fit)
#' }

