% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ridge_learner_no_int.R
\name{makeRidgeLearner_no_int}
\alias{makeRidgeLearner_no_int}
\title{Create a Ridge learner object}
\usage{
makeRidgeLearner_no_int(
  X,
  X_test = NULL,
  growMode = c("NA", "+*", "+", "*"),
  changeToConstant = FALSE
)
}
\arguments{
\item{X}{is our design matrix used in training. Scaling is NOT handled by this function, and must be done by the user beforehand}

\item{X_test}{is our design matrix used in testing Scaling is NOT handled by this function, and must be done by the user beforehand}

\item{growMode}{is a string for if the learner should be grown (or not)
If \code{"+*"}, we grow mu_0 -> (mu_0 * mu_2) + mu_1
If \code{"+"}, we grow mu_0 -> (mu_0 + mu_1)
If \code{"*"}, we grow mu_0 -> (mu_0 * mu_1) (NOTE: Not recommended if we start with \code{k = 1})
If \code{"NA"}, we do not grow this learner}

\item{changeToConstant}{is a logical for if constant fits should be changed to be constant}
}
\description{
Creates a Ridge learner object to be used in \code{\link{veb_boost}}
}
\details{
A Ridge learner \eqn{\beta} has a prior distribution \deqn{\beta ~ N(0, \sigma_{\beta}^2 S^2)} for a fixed and given \eqn{S}.
}
