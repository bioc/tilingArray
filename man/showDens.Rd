\name{showDens}
\alias{showDens}
\title{Plot function for more than one density}
\description{
  Plot function for more than one density
}
\usage{
showDens(z, breaks, col, ylab = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{List of numeric vectors for computing histograms for}
  \item{breaks}{breaks of the histogram}
  \item{col}{colour of the histograms}
  \item{ylab}{y-axis label}
  \item{\dots}{futher arguments passed on to plot}
}
\details{
  ...
}
\value{
  returns scale factor
}
\author{Wolfgang Huber \email{huber@ebi.ac.uk}}
\seealso{\code{\link{hist}}}
\examples{
 # generate 3 vectors of normal-distributed random numbers:
 testdat <- list(x1=rnorm(100), x2=rnorm(100,1,1), x3=rnorm(100,0.5,1))
 # and compare their distribution:
 showDens(testdat, breaks=seq(-4,4,0.5),col=3:5, xlab="Random Numbers")
}
\keyword{internal}% at least one, from doc/KEYWORDS
