\name{plotPenLL}
\alias{plotPenLL}

\title{
  Plot the log-likelihood and penalized
  log-likelihoods (AIC, BIC)}
\description{
  Plot the log-likelihood and two versions of penalized
  log-likelihoods (AIC, BIC) for a segmentation object.
}
  
\usage{
plotPenLL(seg, extrabar=numeric(0), type="b", lty=1, pch=16, lwd=2, ...)
}

\arguments{  
  \item{seg}{A \code{\linkS4class{segmentation}} object.}
  \item{extrabar}{In addition to the location of maximmal BIC,
    vertical bars are drawn at these x-positions as well.}
  \item{type, pch, lty, lwd, ...}{Get passed on to \code{\link{matplot}}.}
}

\details{
  This function is used in the vignette:
  \emph{How to use the segment function to fit a piecewise constant curve.}
}

\value{
  The function is called for its side effect, which is creating a plot
  in the current graphics device.}

\author{Wolfgang Huber <huber@ebi.ac.uk>}

\examples{
  x = rep( sin((0:4)/2*pi), each=3) + rnorm(3*5, sd=0.1)
  res = segment(x, maxseg=8, maxk=15)
  plotPenLL(res)
}
\keyword{manip}
