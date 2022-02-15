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

  if (class(alpha_combined) == "try-error") {
    stop("An error occurred. Be sure that `veb_fit` is the output from a call to `veb_boost_stumps`")
  }

  kl_div_leaves = sapply(veb_fit$leaves, function(x) x$KL_div)
  alpha_combined = lapply(alpha_combined, function(x) if (class(x)[1] == 'try-error') {return(NULL)} else {return(x)})
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
