\name{isWithinInterval}
\alias{isWithinInterval}
\title{Test whether a vector of numbers is within a vector of intervals}
\description{
  Test whether a vector of numbers is within a vector of intervals
}
\usage{
isWithinInterval(x, start, end)
}
\arguments{
  \item{x}{Numeric vector.}
  \item{start}{Numeric vector.}
  \item{end}{Numeric vector of same length as start.}
}
\details{
  For each value \code{x[i]}, this function tests whether it lies
  within the interval from \code{start[j]} to \code{end[j]},
  i.e. whether \code{start[j] <= x[i] <= end[j]}. The result is
  returned in a logical matrix whose rows correspond to \code{i}
  and columns to \code{j}.
  
  The implementation of this function is naive. It does not use loops
  within \code{R}, but can have a large
  memory consumption, corresponding to several matrices of size
  \code{length(x)*length(start)}
}

\value{Matrix of size \code{length(x)*length(start)}.
}

\author{W. Huber <huber@ebi.ac.uk>}

\examples{
  x = 1:8
  start = seq(2, 6, by=2)
  end   = seq(4, 8, by=2)
  isWithinInterval(x, start, end)
}
\keyword{manip}
