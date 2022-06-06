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
#' @param learners is either a single "learner" object, or a list of k "learner" objects
#' A learner object is comprised of:
#' 1. a fit function $fitFunction: (X, Y, sigma2, currentFit) -> newFit (where a fit is a list that must contain $mu1, $mu2, and $KL_div)
#' 2. a prediction function $predFunction: (X, fit, moment) -> posterior moment (1 or 2)
#' 3. a constant check function $constCheckFunction: (fit) -> (TRUE/FALSE) to check if a fit is essentially constant
#' 4. a current fit $currentFit: must contain $mu1 (first posterior moments), $mu2 (second posterior moments), and $KL_div (KL-divergence from q to prior) (can be NULL, at least to start)
#' 5. a predictor object $X (whatever the $fitFunction and $predFunction take in), used for training (can be NULL, e.g. if using constLearner)
#' 6. a predictor object $X_test (whatever the $fitFunction and $predFunction take in), used for testing (can be NULL)
#' 7. a string $growMode for if the learner should be grown (or not)
#' If \code{"+*"}, we grow mu_0 -> (mu_0 * mu_2) + mu_1
#' If \code{"+"}, we grow mu_0 -> (mu_0 + mu_1)
#' If \code{"*"}, we grow mu_0 -> (mu_0 * mu_1) (NOTE: Not recommended if we start with \code{k = 1})
#' If \code{"NA"}, we do not grow this learner
#' 8. a logical $changeToConstant for if the learner should be changed to constant if constCheckFunction evaluates to TRUE
#'
#' @param Y is a numeric response. For all but the 'aft.loglogistic' family, this should be an n-vector.
#' For the 'aft.loglogistic' family, this should be an n x 2 matrix, with the first column being the left end-point of survival time
#' and the second column being the right end-point of survival time (used for interval-censored observations)..
#' In the case of left-censored data, the left end-point should be 'NA', and in the case of right-censored data, the right end-point should be 'NA'.
#' If the observation is uncensored, both end-points should be equal to the observed survival time.
#'
#' @param k is an integer, or a vector of integers of length \code{length(learners)}, for how many terms are in the sum of nodes (for each learner)
#'
#' @param d is either an integer, or an integer vector of length \code{k}, or a list of integer vectors of length \code{length(learners)}
#' (each element either an integer, or a vector of length \code{k}) for the multiplicative depth of each of the k terms
#' NOTE: This can be dangerous. For example, if the fit starts out too large, then entire branhces will be fit to be exactly
#' zero. When this happens, we end up dividing by 0 in places, and this results in NAs, -Inf, etc. USE AT YOUR OWN RISK
#'
#' @param sigma2 is a scalar/n-vector specifying a fixed residual variance. If not NULL, then the residual variance will be
#' fixed to this value/vector. If NULL, then it will be initialized and updated automatically.
#' This should be left as NULL unless you really know what you're doing. For safety, this can only be not NULL if family is gaussian
#'
#' @param family is what family the response is
#'
#' @param tol is a positive scalar specifying the level of convergence to be used
#'
#' @param weights is a vector of the same length as Y weighting the observations in the log-likelihood
#'
#' @param scaleWeights is a logical for if the weights should be scaled to have mean 1 (recommended). If you choose to not scale the weights, then
#' the relative importance of the KL-divergence term will change (possibly desirable, i.e. increase or decrease shrinkage towards the prior)
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
#'
#' learner = list(fitFunction = fitFn, predFunction = predFn, constCheckFunction = constCheckFn, currentFit = NULL, X = X, X_test = NULL, growMode = "+*)
#' veb.fit = veb_boost(learner, family = "gaussian")
#'
#' @export
#'

veb_boost = function(learners, Y, k = 1, d = 1, sigma2 = NULL,
                     family = c("gaussian", "binomial", "multinomial.bouchard", "multinomial.titsias", "negative.binomial", "poisson.log1pexp", "aft.loglogistic", "ordinal.logistic"),
                     weights = 1, scaleWeights = TRUE, exposure = NULL, tol = NROW(Y) / 10000, verbose = TRUE, maxit = Inf, backfit = FALSE) {

  ### Check Inputs ###
  # Family
  family = match.arg(family)
  if (grepl("multinomial", family)) {
    Y = as.character(Y)
  }
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
  # Scalars
  if (any(k < 1) || any(k %% 1 != 0)) {
    stop("'k' must be a positive whole number")
  }
  if (!(length(k) %in% c(1, length(learners)))) {
    stop("'k' must be of length 1 or 'length(learners)'")
  }
  if (any(unlist(d) < 1) || any(unlist(d) %% 1 != 0)) {
    stop("'d' must be a positive whole number")
  }
  # Logical Flags
  if (!(verbose %in% c(TRUE, FALSE))) {
    stop("'verbose' must be either TRUE or FALSE")
  }
  if (!(backfit %in% c(TRUE, FALSE))) {
    stop("'backfit' must be either TRUE or FALSE")
  }
  if (!(scaleWeights %in% c(TRUE, FALSE))) {
    stop("'scaleWeights' must be either TRUE or FALSE")
  }
  # learner object
  if (class(learners) != "list") {
    stop("'learners' must be a list of learners")
  }
  for (learner in learners) {
    if (class(learner$fitFunction) != "function") {
      stop("'$fitFunction' must be a function")
    }
    if (!identical(tolower(names(formals(learner$fitFunction))), c("x", "y", "sigma2", "currentfit"))) {
      stop("'$fitFunction' must take in arguments 'X', 'Y', 'sigma2', and 'currentfit' (in that order, case insensitive)")
    }
    if (class(learner$predFunction) != "function") {
      stop("'$predFunction' must be a function")
    }
    if (!identical(tolower(names(formals(learner$predFunction))), c("x_test", "currentfit", "moment"))) {
      stop("'$predFunction' must take in arguments 'X_test', 'currentFit', and 'moment' (in that order, case insensitive)")
    }
    if (class(learner$constCheckFunction) != "function") {
      stop("'$constCheckFunction' must be a function")
    }
    if (tolower(names(formals(learner$constCheckFunction))) != "currentfit") {
      stop("'constCheckFunction' must take in only the argument currentFit (case insensitive)")
    }
    if (!is.null(learner$growMode) && !(learner$growMode %in% c("+*", "+", "*", "NA"))) {
      stop("'$growMode' must be either NULL, or in c('+*', '+', '*', 'NA')")
    }
  }
  # weights
  if (any(weights <= 0) || any(is.na(weights))) {
    stop("All weights must be positive and not missing")
  }
  if (!(length(weights) %in% c(1, nrow(as.matrix(Y))))) {
    stop("Weights must be of the same length as Y")
  }
  if (length(weights) == 1) {
    weights = rep(weights, nrow(as.matrix(Y)))
  }
  if (scaleWeights) { # scale so the mean is 1
    weights = weights / mean(weights)
  }
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
  if (grepl("multinomial", family)) {
    classes = sort(unique(Y))
    attr(Y, 'which') = sapply(Y, function(j) which(j == classes))
  }

  ### Run Non-multinomial Cases ###
  # if (family != "multinomial") {
  if (!grepl("multinomial", family)) {
    # initialize tree
    veb_boost_tree = initialize_veb_boost_tree(learners = learners, Y = Y, k = k, d = d, weights = weights, family = family, exposure = exposure)
    if (family == "gaussian") {
      veb_boost_tree$sigma2 = if (is.null(sigma2)) var(Y) else sigma2
    } else if (family == "aft.loglogistic") {
      veb_boost_tree$sigma2 = var(rowMeans(attr(Y, 'log_Y'), na.rm = TRUE))
    } else if (family == "ordinal.logistic") {
      veb_boost_tree$cutpoints = seq(from = -max(Y), to = max(Y), length.out = max(Y) - 1)
      attr(veb_boost_tree$cutpoints, "log_Y") = matrix(NA, nrow = length(Y), ncol = 2)
      attr(veb_boost_tree$cutpoints, "log_Y")[attr(Y, 'left_cens'), 2] = veb_boost_tree$cutpoints[1]
      for (k in 2:(length(veb_boost_tree$cutpoints))) {
        attr(veb_boost_tree$cutpoints, "log_Y")[which(Y == k), 1] = veb_boost_tree$cutpoints[k-1]
        attr(veb_boost_tree$cutpoints, "log_Y")[which(Y == k), 2] = veb_boost_tree$cutpoints[k]
      }
      attr(veb_boost_tree$cutpoints, "log_Y")[attr(Y, 'right_cens'), 1] = tail(veb_boost_tree$cutpoints, 1)
    }
    update_sigma2 = ifelse(is.null(sigma2), family %in% c("gaussian", "aft.loglogistic", "ordinal.logistic"), FALSE) # only update variance if Gaussian , or AFT model

    # converge fit once
    veb_boost_tree$convergeFit(tol, update_sigma2 = update_sigma2, verbose = FALSE, maxit = maxit)

    # if growing tree, continue w/ growing tree & fitting to convergence
    if (verbose) {
      cat(paste("ELBO: ", round(veb_boost_tree$ELBO, 3), sep = ""))
      cat("\n")
    }

    veb_boost_tree$convergeFitAll(tol = tol, update_sigma2 = update_sigma2, verbose = FALSE, maxit = maxit)

    # while ((tail(tail(veb_boost_tree$ELBO_progress, 1)[[1]], 1) - tail(tail(veb_boost_tree$ELBO_progress, 2)[[1]], 1) > tol) &&
    #        (length(Traverse(veb_boost_tree, filterFun = function(node) node$isLeaf & !node$isLocked)) > 0) && (tail(tail(veb_boost_tree$ELBO_progress, 1)[[1]], 1) != Inf)) {
    while ((mean(diff(sapply(tail(veb_boost_tree$ELBO_progress, 3), tail, 1))) > 10*tol) &&
           (length(Traverse(veb_boost_tree, filterFun = function(node) node$isLeaf & !node$isLocked)) > 0) && (tail(tail(veb_boost_tree$ELBO_progress, 1)[[1]], 1) != Inf)) {
      if (verbose) {
        cat(paste("ELBO: ", round(veb_boost_tree$ELBO, 3), sep = ""))
        cat("\n")
      }
      veb_boost_tree$convergeFitAll(tol = tol, update_sigma2 = update_sigma2, verbose = FALSE, maxit = maxit)
    }

    veb_boost_tree$lockLearners() # change to constant at the end if needed
    if (backfit) { # converge again, if backfitting
      veb_boost_tree$convergeFit(tol, update_sigma2 = update_sigma2, verbose = FALSE, maxit = Inf)
    }

    veb_boost_tree$predict(1)

    return(veb_boost_tree)

  } else { # else, multinomial case
    classes = sort(unique(Y))
    learner_multiclass = VEBBoostMultiClassLearner$new()
    learner_multiclass$family = family
    learner_multiclass$Y = Y
    learner_multiclass$classes = classes
    fam = ifelse(family == "multinomial.bouchard", "binomial", "multinomial.titsias")

    learnerList = list()
    for (j in 1:length(classes)) { # add learners to learnerList
      classLearner = initialize_veb_boost_tree(learners = learners, Y = 1 * (Y == classes[j]), k = k, d = d, weights = weights, family = fam, my_class_index = j)
      learnerList = c(learnerList, classLearner)
    }
    names(learnerList) = paste0("classLearner_", classes, sep = "")
    learner_multiclass$addclassLearners(learnerList) # add learner list to multi-class clearner

    # converge fit once
    learner_multiclass$convergeFit(tol = tol)

    # if growing tree, continue w/ growing tree & fitting to convergence
    if (verbose) {
      cat(paste("ELBO: ", round(learner_multiclass$ELBO, 3), sep = ""))
      cat("\n")
    }

    learner_multiclass = learner_multiclass$convergeFitAll(tol = tol)

    while ((abs(tail(tail(learner_multiclass$ELBO_progress, 1)[[1]], 1) - tail(tail(learner_multiclass$ELBO_progress, 2)[[1]], 1)) > tol) &&
           (sum(sapply(learner_multiclass$classLearners, function(x) length(Traverse(x, filterFun = function(node) node$isLeaf & !node$isLocked)))) > 0)) {
      if (verbose) {
        cat(paste("ELBO: ", round(learner_multiclass$ELBO, 3), sep = ""))
        cat("\n")
      }
      learner_multiclass$convergeFitAll(tol = tol)
    }

    learner_multiclass$predict(1)

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
#' @param learners is a list of other learners to be used in \code{\link{veb_boost}}
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
#' @param k is an integer, or a vector of integers of length \code{length(learners)}, for how many terms are in the sum of nodes (for each learner)
#'
#' @param d is either an integer, or an integer vector of length \code{k}, or a list of integer vectors of length \code{length(learners)}
#' (each element either an integer, or a vector of length \code{k}) for the multiplicative depth of each of the k terms
#' NOTE: This can be dangerous. For example, if the fit starts out too large, then entire branhces will be fit to be exactly
#' zero. When this happens, we end up dividing by 0 in places, and this results in NAs, -Inf, etc. USE AT YOUR OWN RISK
#'
#' @param use_quants is a logical for if the cut-points should be based off of the quantiles (`use_quants = TRUE`), or if the
#' cut points should be evenly spaced in the range of the variable (`use_quants = FALSE`).
#'
#' @param scale_X is a string for if/how the columns of X should be scaled.
#' 'sd' scales by the standard deviations of the variables.
#' 'max' scales by the maximum absolute value (so variables are on the [-1, +1] scale).
#' 'NA' performs no scaling.
#'
#' @param growMode is a string for if the learner should be grown (or not)
#' If \code{"+*"}, we grow mu_0 -> (mu_0 * mu_2) + mu_1
#' If \code{"+"}, we grow mu_0 -> (mu_0 + mu_1)
#' If \code{"*"}, we grow mu_0 -> (mu_0 * mu_1) (NOTE: Not recommended if we start with \code{k = 1})
#' If \code{"NA"}, we do not grow this learner
#'
#' @param changeToConstant is a logical for if constant fits should be changed to be constant
#'
#' @param max_log_prior_var is a scalar for the maximum that the estimated log-prior variance for each weak learner can be.
#' The idea is that setting this to be small limits the "size" of each weak learner, similar-ish to the learning rate in boosting.
#' The maximum allowed value is 35, which essentially allows the unrestricted MLE to be estimated.
#' The minimum allowed value is -5.
#'
#' @param use_optim is a logical. If TRUE, then the prior variance is optimized using the Brent method.
#' If FALSE, then a single EM step is taken to optimize over V
#'
#' @param lin_prior_prob is a number between 0 and 1 that gives the prior probability that the effect variable is a linear term.
#' This means that (1 - lin_prior_prob) is the probability that the effect variable is a stump term.
#' Within linear terms and stump terms, all variables have the same prior variance.
#'
#' @param reverse_learners is a logical for if the order of learners should be reversed. If FALSE, all additional learners in `learners` will
#' will come first. If TRUE, the stumps learners will be first.
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
#'
#' @useDynLib VEB.Boost, .registration=TRUE
#'
#' @export
#'

veb_boost_stumps = function(X, Y, X_test = NULL, learners = NULL, include_linear = NULL, include_stumps = NULL,
                            num_cuts = ceiling(min(NROW(Y) / 5, max(100, sqrt(NROW(Y))))), k = 1, d = 1,
                            use_quants = TRUE, scale_X = c("sd", "max", "NA"), growMode = c("+*", "+", "*", "NA"), changeToConstant = FALSE,
                            max_log_prior_var = 0, use_optim = TRUE, lin_prior_prob = 0.5, reverse_learners = FALSE, nthreads = ceiling(parallel::detectCores(logical = TRUE) / 2), ...) {
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
  if (any(is.na(X))) {
    stop("'X' cannot have any NAs. Please handle them accordingly first")
  }
  if (!is.null(X_test) && any(is.na(X_test))) {
    stop("'X_test' cannot have any NAs. Please handle them accordingly first")
  }
  if (NROW(X) != NROW(Y)) {
    stop("'X' and 'Y' must have the same number of observations/rows")
  }
  p = ncol(X)
  # check other learners
  if (!is.null(learners) && !is.list(learners)) {
    stop("'learners' must be either NULL or a list of other learners to use")
  }
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
  growMode = match.arg(growMode)
  if (!(changeToConstant %in% c(TRUE, FALSE))) {
    stop("'changeToConstant' must be either TRUE or FALSE")
  }
  # Check max_log_prior_var
  if (!is.numeric(max_log_prior_var) || (length(max_log_prior_var) != 1)) {
    stop("'max_log_prior_var' must be a single number")
  }
  if ((max_log_prior_var < -5) | (max_log_prior_var > 35)) {
    warning("Restricting 'max_log_prior_var' to be in the range [-5, 35]")
    max_log_prior_var = min(35, max(-5, max_log_prior_var))
  }
  if (!(use_optim %in% c(TRUE, FALSE))) {
    stop("'use_optim' must be either TRUE or FALSE")
  }
  # Check lin_prior_prob
  if (!is.numeric(lin_prior_prob) || (length(lin_prior_prob) != 1)) {
    stop("'lin_prior_prob' must be a single number")
  }
  if ((lin_prior_prob < 0) | (lin_prior_prob > 1)) {
    warning("Restricting 'lin_prior_prob' to be in the range [0, 1]")
    lin_prior_prob = min(1, max(0, lin_prior_prob))
  }
  if (!(reverse_learners %in% c(TRUE, FALSE))) {
    stop("'reverse_learners' must be either TRUE or FALSE")
  }

  scale_X = ifelse(scale_X == "sd", 1, ifelse(scale_X == "max", 2, 0))

  if (!is.null(include_linear) && length(include_linear) == 1) {
    include_linear = rep(include_linear, p)
  }
  if (!is.null(include_stumps) && length(include_stumps) == 1) {
    include_stumps = rep(include_stumps, p)
  }

  if (is.infinite(num_cuts) | (!is.null(include_stumps) && all(!include_stumps))) {
    num_cuts = 0
  }
  if (inherits(X, "dgCMatrix")) {
    X_stumps = make_stumps_matrix_sp_cpp(X, include_linear, include_stumps, num_cuts, use_quants, scale_X, nthreads)
  } else {
    X_stumps = make_stumps_matrix_cpp(X, include_linear, include_stumps, num_cuts, use_quants, scale_X, nthreads)
  }

  # set up testing data, if any
  X_test_stumps = NULL
  if (!is.null(X_test)) {
    # if (inherits(X_test, "dgCMatrix") != inherits(X, "dgCMatrix")) {
    #   stop("`X` and `X_test` must both either be sparse (i.e. `dgCMatrix`) or dense. Please convert types to match")
    # }
    if (inherits(X_test, "dgCMatrix")) {
      X_test_stumps = make_stumps_test_matrix_sp_cpp(X_test, X_stumps)
    } else {
      X_test_stumps = make_stumps_test_matrix_cpp(X_test, X_stumps)
    }
  }

  fitFnSusieStumps_maxlV = function(X, Y, sigma2, currentFit) {
    return(weighted_SER_cpp(X, Y, sigma2, currentFit, max_log_prior_var, lin_prior_prob, use_optim))
  }

  stumps_learner = list(fitFunction = fitFnSusieStumps_maxlV, predFunction = predFnSusieStumps_cpp, constCheckFunction = constCheckFnSusieStumps,
                        currentFit = NULL, X = X_stumps, X_test = X_test_stumps, growMode = growMode, changeToConstant = changeToConstant)

  learners = c(learners, list(stumps_learner))
  if (reverse_learners) { # have to explicitly use k and d here, otherwise it doesn't get changed in the `...`
    learners = rev(learners)
    k = rev(k)
    d = rev(d)
  }
  # Run
  veb.fit = veb_boost(learners = learners, Y = Y, k = k, d = d, ...)

  return(veb.fit)
}

