\name{otherStrand}
\alias{otherStrand}
\title{Return the name of the opposite strand}
\description{Return the name of the opposite strand}
\usage{
otherStrand(x)
}
\arguments{
  \item{x}{Character vector whose elements are "+" or "-".}
}
\details{
  This is a rather trivial convenience function.

  An alternative would be to code strands with integers \code{-1} and
  \code{+1}, in which case the inversion would be a trivial builtin
  operation. However, many genomic databases and input data files use
  the character string / factor notation.
}

\value{
  Character vector of same length as \code{x}, with strands reversed.
}

\author{W. Huber <huber@ebi.ac.uk>}

\examples{
  otherStrand(c("+", "-"))
}
\keyword{manip}
