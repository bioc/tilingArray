\name{plot.segmentation}
\alias{plot.segmentation}
\title{Plot result of segmentation algorithm}
\description{
  This function visualizes the result of \code{\link{findSegments}}.
  The user can choose whether he wants to specify the number of segments
  or let be determined by the function (see details).
}
\usage{
  \method{plot}{segmentation}(x, nSegments = NULL, bcol = NULL, from=NULL, to=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class \code{"segmentation"}}
  \item{nSegments}{number of segments in the data, has to be specified by the user}
   %if \code{NULL}, it is determined within the function (see details)}
  \item{bcol}{color(s) to use for segment borders; if \code{NULL} they
    are taken to be the colors 2:\code{nSegments}}
  \item{from}{numeric, index of data point to start plot at; default
    \code{NULL} means take the first one}
  \item{to}{numeric, index of last data point to plot; default
    \code{NULL} means take all up to the last one}
  \item{\dots}{further graphical parameters passed on to
    \code{\link[base]{plot.default}}}
}
\details{
  If \code{nSegments} is not specified, it is taken to be that number,
  at which the Residual Sum of Squares shows the largest decrease from
  the previous one. (The RSS decreases steadily, so absolute minimum is
  reached when \code{nSegments} equals the number of observations.
  This can only be done if confidence intervals have been computed for
  at least one of the segmentations, using the function
  \code{\link{confint.segmentation}}.
}
\value{
  If result is assigned, returns a list containing
  \item{breakp}{the breakpoints of the chosen or computed segmentation}
  \item{confInt}{the confidence intervals of those breakpoints}
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link{findSegments}}, \code{\link{confint.segmentation}}}
\examples{
  dat <- c(rnorm(10,0,1),rnorm(20,2,1), rnorm(5,0.5,0.5),
           rnorm(10,1,1), rnorm(20,2,1))

  if (require("strucchange")){
    segres <- findSegments(dat, maxcp=10, maxk=15, verbose=TRUE)
    segres <- confint(segres, 3:6)
    plot(segres)
  } else {
    segres <- findSegments(dat, maxcp=10, maxk=15, verbose=TRUE) 
    plot(segres, 5)
  }
}
\keyword{hplot}


