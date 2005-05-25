\name{sampleStep}
\alias{sampleStep}
\title{Given a vector of ascending numbers and a step width, sample the
  numbers to keep a certain spacing}
\description{
  Given a vector of ascending numbers and a step width, the function
  determines a subset of those numbers such that there are no
  consecutive triplet of numbers is within an interval \code{step}.
}
\usage{
sampleStep(x, step)
}
\arguments{
  \item{x}{Numeric or integer vector.}
  \item{step}{Numeric scalar.}
}
\details{None yet.}
\value{
  A logical vector of the same length as \code{x}, representing the
  selected subsample.
}

\author{W. Huber <huber@ebi.ac.uk>}

\examples{
 x = sort(as.integer(runif(20)*100))
 sel = sampleStep(x, step=10)
 table(diff(x))
 table(diff(x[sel]))

}
\keyword{manip}
