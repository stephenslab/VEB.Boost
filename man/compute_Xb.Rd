% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_multiplication.R
\name{compute_Xb}
\alias{compute_Xb}
\title{Computes standardized.X \%*\% b}
\usage{
compute_Xb(X, b)
}
\arguments{
\item{X}{an n by p matrix with three attributes: scaled:center, scaled:scale, and attr(X, 'd')}

\item{b}{a p vector}
}
\value{
an n vector
}
\description{
Computes standardized.X \%*\% b
}
\keyword{internal}
