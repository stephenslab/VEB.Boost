#' Initialize VEB-Boost Tree structure
#'
#' Initializes a VEB-Boost tree object as the sum of products of nodes,
#' where you can specify how many learners to add, and the multiplicative depth of each learner.
#'
#' @param X is a prediction objects to be used
#'
#' @param Y is a numeric vector response
#'
#' @param k is an intetger for how many terms are in the sum of nodes
#'
#' @param d is either an integer, or an integer vector of length `k` for the multiplicative depth of each of the k terms
#'
#' @param fitFunctions is a list of length 1 or `k` of fitting functions to be used in
#' each term on the sum of nodes
#'
#' @param predFunctions is a list of length 1 or `k` of prediction functions to be used in
#' each term of the sum of nodes
#'
#' @param constCheckFunctions is a list of length 1 or `k` of constant check functions to be used in
#' each term of the sum of nodes
#'
#' @param family is what family the response is

initialize_veb_boost_tree = function(X, Y, k = 1, d = 1, fitFunctions = list(fitFnSusieStumps), predFunctions = list(predFnSusieStumps),
                                     constCheckFunctions = list(constCheckFnSusieStumps), family = c("gaussian", "binomial")) {
  family = match.arg(family)
  if (length(d) == 1) {
    d = rep(d, k)
  }
  if (length(fitFunctions) == 1) {
    fitFunctions = rep(fitFunctions, k)
  }
  if (length(predFunctions) == 1) {
    predFunctions = rep(predFunctions, k)
  }
  if (length(constCheckFunctions) == 1) {
    constCheckFunctions = rep(constCheckFunctions, k)
  }

  # start by making overall addition of k learners structure
  veb_boost_learner = VEBBoostNode$new(paste("mu_", 0, sep = ""), fitFunction = fitFunctions[[1]], predFunction = predFunctions[[1]], constCheckFunction = constCheckFunctions[[1]], currentFit = list(mu1 = 0, mu2 = 0, KL_div = 0))
  # also add family and response
  veb_boost_learner$family = family
  veb_boost_learner$Y = Y
  for (i in 1:ceiling(log2(k))) {
    base_learners = veb_boost_learner$leaves
    for (learner in base_learners) {
      if (learner$root$leafCount >= k) {
        break
      }
      add_learner = VEBBoostNode$new(paste("mu_", learner$root$leafCount, sep = ""), fitFunction = fitFunctions[[learner$root$leafCount + 1]], predFunction = predFunctions[[learner$root$leafCount + 1]], constCheckFunction = constCheckFunctions[[learner$root$leafCount + 1]], currentFit = list(mu1 = 0, mu2 = 0, KL_div = 0))
      learner$AddSiblingVEB(add_learner, "+", paste("combine_", learner$root$leafCount, sep = ""))
    }
  }

  # now, add predictor object to the root
  veb_boost_learner$X = X

  # now, add multiplicative components, where left-most moments are initialized to 0, and others are initialized to 1 (to avoid infinite variance issue)
  base_learners = veb_boost_learner$leaves
  for (branch in base_learners) {
    j = as.numeric(gsub("mu_", "", branch$name)) + 1 # index of inputs to use for this learner
    depth = d[j]
    mult_count = 1
    for (i in 1:ceiling(log2(depth))) {
      branch_learners = branch$leaves
      for (learner in branch_learners) {
        if (mult_count >= depth) { # if we've added all multiplicative nodes for this branch, go to the next branch
          break
        }
        mult_learner = VEBBoostNode$new(paste("mu_", learner$root$leafCount, sep = ""), fitFunction = fitFunctions[[j]], predFunction = predFunctions[[j]], constCheckFunction = constCheckFunctions[[j]], currentFit = list(mu1 = 1, mu2 = 1, KL_div = 0))
        learner$AddSiblingVEB(mult_learner, "*", paste("combine_", learner$root$leafCount, sep = ""))
        mult_count = mult_count + 1
      }
    }
  }

  return(veb_boost_learner)

}
