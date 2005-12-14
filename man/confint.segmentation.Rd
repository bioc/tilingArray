\name{confint.segmentation}
\alias{confint.segmentation}
\title{Confidence intervals for a segmentation fit}
\description{
  Wrapper for \code{\link{confint.breakpointsfull}} from the
  \code{strucchange} package for
  computing asymptotic confidence intervals for a segmentation fit.
}

\usage{
\method{confint}{segmentation}(object, n, ...)
}

\arguments{
  \item{object}{object of class \code{"segmentation"};
    e.g. the result from \code{\link{findSegments}}.}
  \item{n}{numeric vector, selects the segmentation fits
    for which to compute confidence intervals for.
    Typically there are multiple such fits, one for each choice of total
    number of segments, see \code{\link{findSegments}}.}
  \item{...}{passed on to \code{\link{confint.breakpointsfull}}.}
}

\details{
  See \code{\link{confint.breakpointsfull}}
}

\value{
  An object of class \code{"segmentation"}, namely a copy of the argument
  \code{object}, with added slots for the confidence intervals.
}

\author{Wolfgang Huber \email{huber@ebi.ac.uk}}
\seealso{
  \code{\link{segment}},
  \code{\link{plot.segmentation}},
  \code{\link[stats]{confint}},
  \code{\link[strucchange]{confint.breakpointsfull}}
}

\examples{
   # generate random data with 5 segments:
   x <- c(rnorm(10,0,1),rnorm(10,3,1),rnorm(10,0.5,0.5),
          rnorm(10,1.5,1),rnorm(10,5,1))

   segres <-  findSegments(x, maxcp=10, maxk=15)
   segres <- confint.segmentation(segres,nSegments=c(3,4,5,6))

   # see that the step between segments 3 and 4 is less certain than
   #  the other ones:
   segres$confInt
   plot(segres,5, pch=20)
}
\keyword{manip}

