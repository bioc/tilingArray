\name{breakpointsPretend}
\alias{residuals.breakpointsPretend}
\alias{breakpoints.breakpointsPretend}
  
\title{Accessor methods for breakpointsPretend objects - not to be called
  by the user.}

\description{
  Accessor methods for breakpointsPretend objects - not to be called
  by the user.

  These functions are used in the interface between the
  \code{segmentation} class and the
  \code{\link[strucchange:confint.breakpointsfull]{confint.breakpointsfull}}
  method of the strucchange package. This method calls
  \code{breakpoints} and \code{residuals} methods for its first
  argument, and since we pass an argument of S3 class
  \code{breakpointsPretend}, we can avoid the overhead of the
  corresponding methods for \code{breakpointsfull} objects.
  These functions are of no interest to the user.
}

\usage{
\S3method{residuals}{breakpointsPretend}(object, breaks, ...)
breakpoints.breakpointsPretend(obj, breaks, ...) 
}

\arguments{
  \item{object}{a breakpointsPretend object.}
  \item{obj}{a breakpointsPretend object.}
  \item{breaks}{dummy argument, is ignored.}
  \item{...}{futher arguments.}
}

\value{residuals and breakpoints.}

\author{W. Huber \email{huber@ebi.ac.uk}}

\keyword{manip}
