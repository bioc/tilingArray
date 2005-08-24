\name{confint.segmentation}
\alias{confint.segmentation}
\title{Compute Bai's confidence intervals for specified segmentations}
\description{
  Compute Bai's confidence intervals for specified segmentations
}
\usage{
\method{confint}{segmentation}(object,parm="breakpoints", level = 0.95,
                                 nSegments=NULL, 
                                 het.reg = FALSE,het.err = FALSE, ...)

}
\arguments{
  \item{object}{Object of class \code{"segmentation"};
    Result of \code{findSegments}}
  \item{parm}{character; which parameters to compute confidence
    intervals for; only "breakpoints" is implemented}
  \item{level}{Confidence level of the confidence intervals}
  \item{nSegments}{Number of segments in segmentations to compute
    confidence intervals for. Defaults to computing them for all
    segmentations considered in \code{findSegments}}
  \item{het.reg}{logical. Should heterogenous regressors be assumed? If set
    to \code{FALSE} the distribution of the regressors is assumed to be
    homogenous over the segments.}
  \item{het.err}{logical. Should heterogenous errors be assumed? If set
    to \code{FALSE} the distribution of the errors is assumed to be
    homogenous over the segments.}
  \item{...}{currently not used}
}
\details{
  Basically, this function just prepares an object for calling
  the function \code{\link{computeConfInt}}.
  
  The distribution function used for the computation of confidence
  intervals of breakpoints is given in Bai (1997). The procedure, in
  particular the usage of heterogenous regressors and/or errors, is
  described in more detail in Bai & Perron (2003).

  The breakpoints should be computed from a formula with \code{breakpoints},
  then the confidence intervals for the breakpoints can be derived by
  \code{confint} and these can be visualized with the segmentation.
  For an example see \code{\link{plot.segmentation}}.

}
\value{
  An object of class \code{"segmentation"}. Actually the same as the
  argument \code{object} with the following list items replaced
  \item{chosenSegNo}{Numeric; Segment numbers of segmentations,
    for which confidence intervals were computed}
  \item{confInt}{List of confidence intervals as tables for those
    segmentations}
  \item{residuals}{List of numeric vectors containing the residuals for
    those segmentations}
  \item{call}{with call of this function appended}
}
\references{
  Bai J. (1997), Estimation of a Change Point in Multiple Regression Models,
  \emph{Review of Economics and Statistics}, \bold{79}, 551-563.

  Bai J., Perron P. (2003), Computation and Analysis of Multiple
  Structural Change Models,
  \emph{Journal of Applied Econometrics}, \bold{18}, 1-22.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link{computeConfInt}},\code{\link{findSegments}},
  \code{\link{plot.segmentation}},\code{\link[stats]{confint}}}
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

