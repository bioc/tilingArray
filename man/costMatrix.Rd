\name{costMatrix}
\alias{costMatrix}
\title{Segmentation cost matrix}
\description{
  This function calculates the cost matrix for the segmentation model
}
\usage{
costMatrix(x, maxk)
}
\arguments{
  \item{x}{Numeric vector of length \code{n} or matrix with \code{n} rows and \code{d}
    columns, where \code{n} is the number of sample points and \code{d}
    the number of replicate measurements (e.g. from multiple arrays).}
  \item{maxk}{Positive integer.}
}
\details{
}

\value{
  Matrix with \code{maxk} rows and \code{length(x)} columns.
}

\author{W. Huber}

\examples{
d = 4
x = apply(matrix(rnorm(200), ncol=d), 2, cumsum)
maxk = 30
 
G = costMatrix(x, maxk=maxk)
 
G.pedestrian = matrix(NA, nrow=nrow(G), ncol=ncol(G))
for(i in 1:(ncol(G)))
  for(k in 1:min(nrow(G), nrow(x)-i+1))
    G.pedestrian[k, i] = (k*d-1)/d*var(as.vector(x[i:(i+k-1), ]))
 
stopifnot(identical(is.na(G), is.na(G.pedestrian)))
stopifnot(max(abs(G-G.pedestrian), na.rm=TRUE) <= 1e-6)
}
\keyword{manip}
