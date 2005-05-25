\name{sampleStep}
\alias{sampleStep}
\title{Sampling of ascending numbers to ensure minimal spacing.}
\description{
  Given a vector of ascending numbers and a step width, sample the
  numbers such that the difference between consecutive numbers is greater
  than or equal to \code{step}.
}
\usage{
sampleStep(x, step)
}
\arguments{
  \item{x}{Numeric or integer vector.}
  \item{step}{Numeric scalar.}
}
\details{The simple algorithm works greedily from \code{x[1]} to
  \code{x[length(x)]}. First, \code{x[1]} is selected. Then, if \code{x[i]}
  is selected, all numbers \code{x[j]} with \code{j>i} and
  \code{x[j]-x[i]<step} are dropped. Then, \code{i} is set to the
  smallest j with \code{x[j]-x[i]>=step}.}
\value{
  A logical vector of the same length as \code{x}, representing the
  selected subsample.
}

\author{W. Huber <huber@ebi.ac.uk>}

\examples{
 x   = sort(as.integer(runif(20)*100))
 sel = sampleStep(x, step=10)
 x
 x[sel]

}
\keyword{manip}
