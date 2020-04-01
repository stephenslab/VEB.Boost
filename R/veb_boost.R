#' Performs VEB-Boosting
#'
#' Solves the VEB-Boost regression problem using the supplied inputs
#'
#' @details
#'
#' Given a pre-specified arithmetic tree structure \deqn{T(\mu_1, \dots, \mu_L)},
#' priors \deqn{\mu_l \sim g_l(\cdot)}, and inputs for the response, VEB-Boosting is performed.
#'
#' A cyclic CAVI scheme is used, where we cycle over the leaf nodes and update the approxiomation
#' to the posterior distribution at each node in turn.
#'
#' We start with the arithmetic tree structure \deqn{T(\mu_1, \dots, \mu_L) = \sum_{i=1}^k \prod_{j=1}^{d_k} \mu_{i, j}}
#'
#'
#' @param X is a list of predictor objects with either 1 or k elements. If it contains 1 element, then the same
#' predictor is used for each base learner. Otherwise, the k sub-trees that we add together each get their own X.
#' Each predictor can take any form, so long as the user-supplied fitFunctions and predFunctions know how to use them.
#'
#' @param Y is a numeric vector response
#'
#' @param fitFunctions is either a single fitting function, or a list of length `k` of fitting functions to be used in
#' each term on the sum of nodes
#'
#' @param predFunctions is either a single prediction function, or a list of length `k` of prediction functions to be used in
#' each term of the sum of nodes
#'
#' @param constCheckFunctions is either a single constant check function, or a list of length `k` of constant check functions
#' to be used in each term of the sum of nodes
#'
#' @param growTree is a logical for if we should grow the tree after convergence (TRUE), or only use the initial tree
#' structure (FALSE)
#'
#' @param k is an integer for how many terms are in the sum of nodes
#'
#' @param d is either an integer, or an integer vector of length `k` for the multiplicative depth of each of the k terms
#' NOTE: This can be dangerous. For example, if the fit starts out too large, then entire branhces will be fit to be exactly
#' zero. When this happens, we end up dividing by 0 in places, and this results in NAs, -Inf, etc. USE AT YOUR OWN RISK
#'
#' @param growMode specifies how we grow the tree, either splitting nodes with addition, multiplication, or both
#' If `+`, we grow mu_0 -> (mu_0 + mu_1)
#' If `*`, we grow mu_0 -> (mu_0 * mu_1) (NOTE: Not recommended if we start with `k = 1`)
#' If `+*``, we grow mu_0 -> (mu_0 * mu_2) + mu_1
#'
#' @param changeToConstant is a flag for if, when the fit is found to be basically constant, if we should actually change
#' the fitting function of that node to fit exactly a constant value
#'
#' @param family is what family the response is
#'
#' @param tol is a positive scalar specifying the level of convergence to be used
#'
#' @param verbose is a logical flag specifying whether we should report convergence information as we go
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
#' predFn = function(X, fit, moment = c(1, 2)) {
#'   # For input X and list `fit` returned from fitFn (encoding the variational posterior for b), computes either the 1st or 2nd 
#'   # posterior moments of the response (depending on if `moment` is 1 or 2)
#'   if (moment == 1) {
#'     res = E_fit[Xb]
#'   } else {
#'     res = E_fit[(Xb)^2]
#'   }
#'   return(res)
#' }
#' fitFn = function(X, Y, sigma2, init) {
#'   # For a given prior g(b), a function that approximates the posterior of b, q(b), using Variational Inference
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
#' constCheckFn = function(fit) {
#'   # For a given `fit`, returns TRUE if the variational posterior is close enough to a constant, else FALSE
#'   if (fit is close to constant) {
#'     return(TRUE)
#'   } else {
#'     return(FALSE)
#'   }
#' }
#' veb.fit = veb_boost(list(X), Y, fitFn, predFn, constCheckFn, family = "gaussian")
#'
#' @export
#'

veb_boost = function(X, Y, fitFunctions, predFunctions, constCheckFunctions,
                     growTree = TRUE, k = 1, d = 1, growMode = c("+", "*", "+*"), changeToConstant = TRUE,
                     family = c("gaussian", "binomial", "multinomial"),
                     tol = length(Y) / 10000, verbose = TRUE, mc.cores = 1) {

  ### Check Inputs ###
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
  # Scalars
  if ((k < 1) || (k %% 1 != 0)) {
    stop("'k' must be a positive whole number")
  }
  if ((d < 1) || (d %% 1 != 0)) {
    stop("'d' must be a positive whole number")
  }
  # Predictors
  if ((class(X) != "list") || !(length(X) %in% c(1, k))) {
    stop("'X' must be a list of length 1 or k")
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
  # mc.cores
  if (mc.cores != 1) {
    message("Parallel multinomial updates not currently supported, setting 'mc.cores' to 1")
    mc.cores = 1
  }

  ### Run Gaussian/Binomial Cases ###
  if (family %in% c("gaussian", "binomial")) {
    # initialize tree
    mu = initialize_veb_boost_tree(Xs = X, Y = Y, k = k, d = d, fitFunctions = fitFunctions, predFunctions = predFunctions,
                                               constCheckFunctions = constCheckFunctions, family = family)
    if (family == "gaussian") {
      mu$sigma2 = var(Y)
    }
    update_sigma2 = (family == "gaussian") # only update variance if Gaussian response

    # converge fit once
    mu$convergeFit(tol, update_sigma2 = update_sigma2, verbose = FALSE)

    # if growing tree, continue w/ growing tree & fitting to convergence
    if (growTree) {
      learner = mu$convergeFitAll(tol = tol, update_sigma2 = update_sigma2, growMode = growMode, changeToConstant = changeToConstant, verbose = FALSE)

      while ((abs(tail(tail(learner$ELBO_progress, 1)[[1]], 1) - tail(tail(learner$ELBO_progress, 2)[[1]], 1)) > tol) && (length(Traverse(learner, filterFun = function(node) node$isLeaf & !node$isLocked)) > 0)) {
        if (verbose) {
          cat(paste("ELBO: ", round(learner$ELBO, 3), sep = ""))
          cat("\n")
        }
        learner$convergeFitAll(tol = tol, update_sigma2 = update_sigma2, growMode = growMode, changeToConstant = changeToConstant, verbose = FALSE)
      }

      return(learner)
    } else {
      return(mu)
    }
  } else { # else, multinomial case
    classes = sort(unique(Y))
    learner_multiclass = VEBBoostMultiClassLearner$new()
    learner_multiclass$Y = Y
    learner_multiclass$mc.cores = mc.cores
    learner_multiclass$classes = classes
    if (length(X) == 1) { # if the same predictor object is to be used for all nodes, just store once in the multi-class learner
      learner_multiclass$X = X[[1]]
      X = NULL
    }
    learnerList = list()
    for (j in 1:length(classes)) { # add learners to learnerList
      learner = initialize_veb_boost_tree(Xs = X, Y = 1 * (Y == classes[j]), k = k, d = d, fitFunctions = fitFunctions, predFunctions = predFunctions,
                                constCheckFunctions = constCheckFunctions, family = "binomial")
      learnerList = c(learnerList, learner)
    }
    names(learnerList) = paste0("learner", 0:(length(classes) - 1), sep = "")
    learner_multiclass$addLearners(learnerList) # add learner list to multi-class clearner

    # converge fit once
    learner_multiclass$convergeFit(tol = tol)

    # if growing tree, continue w/ growing tree & fitting to convergence
    if (growTree) {
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

    return(learner_multiclass)
  }

}


#' Wrapper for using SER w/ stumps
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
#' @param include_linear is a logical indicator specifying whether we should include linear terms from X
#' 
#' @param include_stumps is a logical indicator specifying whether we should include stumps terms from X
#' 
#' @param num_cuts is a whole number specifying how many cuts to make when making the stumps terms.
#' We use the quantiles from each predictor when making the stumps splits, using `num-cuts` of them.
#' If `num_cuts = NULL`, then all values of the variables are used as split points.
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
#' @export
#' 

veb_boost_stumps = function(X, Y, include_linear = TRUE, include_stumps = TRUE, num_cuts = 100, ...) {
  ### Check Inputs ###
  # Check logicals
  if (!(include_linear %in% c(TRUE, FALSE))) {
    stop("'include_linear' must be either TRUE or FALSE")
  }
  if (!(include_stumps %in% c(TRUE, FALSE))) {
    stop("'include_stumps' must be either TRUE or FALSE")
  }
  # Make sure at least one include is true
  if (!(include_linear | include_stumps)) {
    stop("At least one of 'include_linear' or 'include_stumps' must be TRUE")
  }
  # Make sure num_cuts is a positive whole number
  if ((num_cuts < 1) || (num_cuts %% 1 != 0)) {
    stop("'num_cuts' must be a positive whole number")
  }
  
  # Make stumps matrix
  if (is.null(num_cuts)) {
    cuts = NULL
  } else {
    cuts = apply(X, MARGIN = 2, function(col) quantile(col, probs = seq(from = 0, to = 1, length.out = num_cuts)))
    cuts = sapply(1:ncol(cuts), function(i) list(cuts[, i]))
  }
  X_stumps = make_stumps_matrix(X, include_linear, include_stumps, cuts)

  # Run
  return(veb_boost(X = list(X_stumps), Y = Y,
                   fitFunctions = fitFnSusieStumps, predFunctions = predFnSusieStumps, constCheckFunctions = constCheckFnSusieStumps,
                   ...))
}
