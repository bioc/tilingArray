\name{posMin}
\alias{posMin}
\title{Find the smallest positive number in a vector}
\description{Find the smallest positive number in a vector}
\usage{
posMin(x, ...)
}
\arguments{
  \item{x}{Numeric vector.}
  \item{...}{Further arguments that get passed on to \code{min}.}
}
\details{
  This is a rather trivial convenience function.
}

\value{Numeric of length 1.
}

\author{W. Huber <huber@ebi.ac.uk>}

\examples{
  x = runif(5)
  posMin(x-0.5)
  posMin(x-2)
}
\keyword{manip}
