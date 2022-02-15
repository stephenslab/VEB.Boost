#' Initialize VEB-Boost Tree structure
#'
#' Initializes a VEB-Boost tree object as the sum of products of nodes,
#' where you can specify how many learners to add, and the multiplicative depth of each learner.
#'
#' @param learners is either a single "learner" object, or a list of k "learner" objects
#' A learner object is comprised of:
#' 1. a fit function $fitFunction: (X, Y, sigma2, currentFit) -> newFit (where a fit is a list that must contain $mu1, $mu2, and $KL_div)
#' 2. a prediction function $predFunction: (X, fit, moment) -> posterior moment (1 or 2)
#' 3. a constant check function $constCheckFunction: (fit) -> (TRUE/FALSE) to check if a fit is essentially constant
#' 4. a current fit $currentFit: must contain $mu1 (first posterior moments), $mu2 (second posterior moments), and $KL_div (KL-divergence from q to prior) (can be NULL, at least to start)
#' 5. a predictor object $X (whatever the $fitFunction and $predFunction take in), used for training (can be NULL, e.g. if using constLearner)
#' 6. a predictor object $X_test (whatever the $fitFunction and $predFunction take in), used for testing (can be NULL)
#'
#' @param Y is a numeric vector response
#'
#' @param k is an intetger for how many terms are in the sum of nodes
#'
#' @param d is either an integer, or an integer vector of length `k` for the multiplicative depth of each of the k terms
#'
#' @param weights is a vector of the same length as Y weighting the observations. Relative weights are used and we take care of the scaling for you
#'
#' @param family is what family the response is

initialize_veb_boost_tree = function(learners, Y, k = 1, d = 1, weights = 1,
                                     family = c("gaussian", "binomial", "negative.binomial", "poisson.log1pexp", "aft.loglogistic", "ordinal.logistic"), exposure = NULL) {
  family = match.arg(family)
  if (length(d) == 1) {
    d = rep(d, k)
  }
  if (length(learners) == 1) {
    learners = rep(learners, k)
  }

  # 1/k of approximate average
  if (family == "gaussian") {
    mu_init = mean(Y) / k
  } else if (family == "binomial") {
    mu_init = log(mean(Y) / (1 - mean(Y))) / k
  } else if (family == "negative.binomial") {
    mu_init = log(mean(Y / exposure)) / k
  } else if (family == "poisson.log1pexp") {
    mu_init = mean(Y / exposure) / k
    # mu_init = 0
  } else if (family == "poisson.exp") {
    mu_init = mean(log(Y + 1e-8) - log(exposure)) / k
    # mu_init = 0
  } else if (family == "aft.loglogistic") {
    mu_init = mean(log(Y), na.rm = TRUE) / k
  } else {
    mu_init = 0
  }

  # start by making overall addition of k learners structure
  veb_boost_learner = VEBBoostNode$new("mu_0", learner = learners[[1]])
  if (is.null(veb_boost_learner$learner$currentFit)) {
    veb_boost_learner$learner$currentFit = list(mu1 = mu_init, mu2 = mu_init^2, KL_div = 0)
  }
  # also add family and response
  veb_boost_learner$family = family
  veb_boost_learner$Y = Y
  veb_boost_learner$weights = weights
  for (i in 1:ceiling(log2(k))) {
    base_learners = veb_boost_learner$leaves
    for (base_learner in base_learners) {
      if (base_learner$root$leafCount >= k) {
        break
      }
      add_learner = VEBBoostNode$new(paste("mu_", base_learner$root$leafCount, sep = ""), learner = learners[[base_learner$root$leafCount + 1]])
      if (is.null(add_learner$learner$currentFit)) {
        add_learner$learner$currentFit = list(mu1 = mu_init, mu2 = mu_init^2, KL_div = 0)
      }
      base_learner$AddSiblingVEB(add_learner, "+", paste("combine_", base_learner$root$leafCount, sep = ""))
    }
  }

  # and add exposure
  veb_boost_learner$exposure = exposure

  # now, add multiplicative components, where left-most moments are initialized to 0, and others are initialized to 1 (to avoid infinite variance issue)
  base_learners = veb_boost_learner$leaves
  for (branch in base_learners) {
    # if (branch$isLocked) {
    #   next
    # }
    j = as.numeric(gsub("mu_", "", branch$name)) + 1 # index of inputs to use for this learner
    depth = d[j]
    mult_count = 1
    for (i in 1:ceiling(log2(depth))) {
      branch_learners = branch$leaves
      for (branch_learner in branch_learners) {
        if (mult_count >= depth) { # if we've added all multiplicative nodes for this branch, go to the next branch
          break
        }
        mult_learner = VEBBoostNode$new(paste("mu_", branch_learner$root$leafCount, sep = ""), learner = learners[[j]])
        mult_learner$learner$currentFit = list(mu1 = 1, mu2 = 1, KL_div = 0)
        branch_learner$AddSiblingVEB(mult_learner, "*", paste("combine_", branch_learner$root$leafCount, sep = ""))
        mult_count = mult_count + 1
      }
    }
  }

  return(veb_boost_learner)

}
