\name{Gmean}
\alias{Gmean}
\title{Segmentation cost matrix}
\description{
  This function calculates the cost matrix for a homoscedastic model
}
\usage{
Gmean(x, maxk)
}
\arguments{
  \item{x}{Numeric (real) vector.}
  \item{maxk}{Positive integer.}
}
\details{
}

\value{
  Matrix with \code{maxk} rows and \code{length(x)} columns.
}

\author{W. Huber}

\examples{
x = cumsum(rnorm(100))
G = Gmean(x, maxk=20)


G.pedestrian = matrix(NA, nrow=nrow(G), ncol=ncol(G))
G.pedestrian[1,] = 0
for(i in 1:(ncol(G)-1))
  for(k in 2:min(nrow(G), length(x)-i+1))
    G.pedestrian[k, i] = (k-1)*var(x[i:(i+k-1)])

stopifnot(identical(is.na(G), is.na(G.pedestrian)))
stopifnot(max(abs(G-G.pedestrian), na.rm=TRUE) <= 1e-6)

}
\keyword{manip}
