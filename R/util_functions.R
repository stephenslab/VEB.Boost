#' @keywords internal
ilogit = function(x) { # inverse logit function
  1 / (1 + exp(-x))
}

# numerically stable log(1 + exp(x))
# see https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf formula 10
#' @keywords internal
log1pexp = function(x) {
  suppressWarnings(rowSums(cbind((exp(x) * (x <= -37)),
                                 (log1p(exp(x)) * (x > -37) * (x <= 18)),
                                 ((x + exp(-x)) * (x > 18) * (x <= 33.3)),
                                 (x * (x > 33.3))), na.rm = T))
}

# numerically stable log(log(1 + exp(x)))
# see https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf formula 10
#' @keywords internal
loglog1pexp = function(x) {
  suppressWarnings(rowSums(cbind((x * (x <= -37)),
                                 (log(log1p(exp(x))) * (x > -37) * (x <= 18)),
                                 (log(x + exp(-x)) * (x > 18) * (x <= 33.3)),
                                 (log(x) * (x > 33.3))), na.rm = T))
}


#' Calculates variable importance measure from a `veb_boost_stumps` fit
#'
#' This function can calculate a few variable importance measures from a `veb_boost_stumps` fit. It
#' combines the posterior probabilities from each base learner (optionally weighted by the KL-divergence
#' of that base learner) in order to give a single numerical value for variable importance for each variable
#' (combinind linear and stumps terms together)
#'
#' @details
#' This function takes in a return object from \code{veb_boost_stumps} and returns a measure of variable importance.
#' It either returns an estimated PIP for if a variable is included in some form (linear or stump term) somewhere in the fit,
#' or the sum of probabilities across all base learners for that variable (combining linear and stump terms)
#'
#' @param veb_fit is the fitted object from \code{veb_boost_stumps}
#'
#' @param method is either "pip" or "sum", depending on how the posterior probabilities from the base learners are to be combined
#' N.B. "pip" isn't very informative if `veb_fit` has many leaves and `scale_by_KL_div` is set to \code{FALSE}
#'
#' @param scale_by_KL_div is a logical for if the KL-divergence of the base learner should be taken into account.
#' Intuitively, if a learner has a learner KL-divergence, it is fitting more signal, and should count for more
#'
#' @return a vector of positive values, either PIPs or sum of probabilities. The order is in the same order as the columns of the
#' original \code{X} matrix that was supplied to \code{veb_boost_stumps}
#'
#' @export
stumpsVariableImportance = function(veb_fit, method = c("pip", "sum"), scale_by_KL_div = TRUE) {
  method = match.arg(method)
  if (!(scale_by_KL_div %in% c(TRUE, FALSE))) {
    stop("'scale_by_KL_div' must be either TRUE or FALSE")
  }

  alpha_combined = try({
    lapply(veb_fit$leaves, function(x) try({getAlphaByVar(x$learner$X, x$learner$currentFit)}, silent = T))
  }, silent = TRUE)

  if (inherits(alpha_combined, "try-error")) {
    stop("An error occurred. Be sure that `veb_fit` is the output from a call to `veb_boost_stumps`")
  }

  kl_div_leaves = sapply(veb_fit$leaves, function(x) x$KL_div)
  alpha_combined = lapply(alpha_combined, function(x) if (inherits(class(x)[1], 'try-error')) {return(NULL)} else {return(x)})
  is_null_alpha = sapply(alpha_combined, is.null)
  kl_div_leaves = kl_div_leaves[!is_null_alpha]
  alpha_combined = do.call(cbind, alpha_combined)
  # alpha_combined = do.call(cbind, lapply(alpha_combined, function(x) if (class(x)[1] == 'try-error') {return(NULL)} else {return(x)}))

  if (scale_by_KL_div) {
    kl_div_leaves = kl_div_leaves / mean(kl_div_leaves)
  } else {
    kl_div_leaves = rep(1, length(kl_div_leaves))
  }
  kl_div_leaves = matrix(kl_div_leaves, nrow = nrow(alpha_combined), ncol = length(kl_div_leaves), byrow = TRUE)

  if (method == "pip") {
    suppressWarnings({
      a = 1 - exp(rowSums(kl_div_leaves * log(1 - alpha_combined)))
    })
    a[is.nan(a)] = 1
  } else if (method == "sum") {
    a = rowSums(kl_div_leaves * alpha_combined)
  } else {
    stop("Unrecognized `method` value")
  }
  return(a)
}


#' Performs a prediction using a `veb_boost` fit and supplied X_test objects
#'
#' This function performs a prediction using a `veb_boost` fit and supplied X_test objects.
#' If the `veb_boost` object was created with `L` different learners, then `X_test_list` must be
#' a list of length `L`, in the order corresponding to the learners
#'
#' @details
#' This function takes in a return object from \code{\link{veb_boost}} and a list of predictor objects and calculates
#' the desired posterior moments.
#'
#' @param veb_fit is the fitted object from \code{\link{veb_boost}}
#'
#' @param X_test_list is a list of predictor objects. This must be of the same length as the list of learners used
#' in fitting `veb_fit`
#'
#' @param moment is the desired posterior moment to calculate. This can be either 1, 2, or c(1, 2)
#'
#' @return a vector of first posterior moments, second posterior moments, or a matrix whose columns are the first and
#' second posterior moments respectively (depending on the supplied `moment`)
#'
#' @export
predict_veb_boost = function(veb_fit, X_test_list, moment = c(1, 2)) {
  if (!identical(moment, 1) && !identical(moment, 2) && !identical(moment, c(1, 2))) {
    stop("'moment' must be either 1, 2, or c(1, 2)")
  }
  if (!inherits(X_test_list, "list")) {
    stop("'X_test_list' must be a list")
  }
  combine_full_nodes = list()
  # first, clear existing X_test values and get a list of all nodes that combine different types of learners
  veb_fit$Do(function(x) {
    x$learner$X_test = NULL
    if (grepl("combine_full", x$name)) {
      combine_full_nodes <<- c(combine_full_nodes, x)
    }
  }, traversal = 'post-order')
  # check the number of predictor objects matches
  if (length(combine_full_nodes) + 1 != length(X_test_list)) {
    stop(paste("'X_test_list' is of length ", length(X_test_list), ", but should be of length ", length(combine_full_nodes) + 1, sep = ""))
  }
  # now, add new X_test values
  if (length(combine_full_nodes) == 0) {
    try({veb_fit$root$learner$X_test = X_test_list[[1]]}, silent = TRUE)
  } else {
    for (node in combine_full_nodes) {
      j = as.numeric(strsplit(node$name, "combine_full_")[[1]][2]) + 1
      try({node$children[[1]]$learner$X_test = X_test_list[[j]]}, silent = TRUE)
      try({node$children[[2]]$learner$X_test = X_test_list[[j + 1]]}, silent = TRUE)
    }
  }
  # now, predict and return
  veb_fit$predict(moment)
  if (identical(moment, 1)) {
    return(veb_fit$pred_mu1)
  } else if (identical(moment, 2)) {
    return(veb_fit$pred_mu2)
  } else {
    return(cbind(veb_fit$pred_mu1, veb_fit$pred_mu2))
  }
}


#' Performs a prediction using a `veb_boost_stumps` fit and supplied X_test objects (for non-stumps learners)
#'
#' This function performs a prediction using a `veb_boost_stumps` fit and supplied X_test objects.
#' If the `veb_boost_stumps` object was created with `L` different learners, then `X_test_list` must be
#' a list of length `L-1`, in the order corresponding to the learners. The last predictor object is made
#' using the stumps matrix from the stumps learners
#'
#' @details
#' This function takes in a return object from \code{\link{veb_boost_stumps}} and a list of predictor objects and calculates
#' the desired posterior moments.
#'
#' @param veb_fit_stumps is the fitted object from \code{\link{veb_boost_stumps}}
#'
#' @param X_test is a matrix to be used when making the new stumps matrix for the stumps learners
#'
#' @param X_test_list is a list of predictor objects. This must be of the same length as the list of learners used
#' in fitting `veb_fit_stumps`
#'
#' @param moment is the desired posterior moment to calculate. This can be either 1, 2, or c(1, 2)
#'
#' @return a vector of first posterior moments, second posterior moments, or a matrix whose columns are the first and
#' second posterior moments respectively (depending on the supplied `moment`)
#'
#' @export
predict_veb_boost_stumps = function(veb_fit_stumps, X_test, X_test_list = NULL, moment = c(1, 2)) {
  if (!is_valid_matrix(X_test)) {
    stop("'X_test' must be a valid numeric matrix")
  }
  if (!is.null(X_test_list) && class(X_test_list) != "list") {
    stop("'X_test_list' must be NULL or a list")
  }
  X_test_stumps = NULL
  for (node in rev(veb_fit_stumps$leaves)) {
    if (body(node$learner$predFunction) == body(predFnSusieStumps_cpp)) {
      if (inherits(X_test, "dgCMatrix")) {
        X_test_stumps = make_stumps_test_matrix_sp_cpp(X_test, node$learner$X)
      } else {
        X_test_stumps = make_stumps_test_matrix_cpp(X_test, node$learner$X)
      }
      break
    }
  }
  if (is.null(X_test_stumps)) {
    stop("No stumps learner found in 'veb_fit_stumps")
  }
  return(predict_veb_boost(veb_fit_stumps, c(X_test_list, X_test_stumps), moment))
}
