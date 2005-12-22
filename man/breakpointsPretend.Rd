\name{breakpointsPretend}
\alias{residuals.breakpointsPretend}
\alias{breakpoints.breakpointsPretend}
  
\title{Accessor methods for breakpointsPretend objects - not to be called
  by the user.}
\description{
  Accessor methods for breakpointsPretend objects - not to be called
  by the user. These functions are needed for the interface between the
  \code{segmentation} class and the
  \code{\link[strucchange:confint.breakpointsfull]{confint.breakpointsfull}}
  method of the strucchange package. They are of no interest to the
  user.
}

\usage{
residuals.breakpointsPretend(object, breaks, ...)
breakpoints.breakpointsPretend(obj, breaks, ...) 
}

\arguments{
  \item{object}{a breakpointsPretend object.}
  \item{obj}{a breakpointsPretend object.}
  \item{...}{futher arguments.}
}

\value{residuals and breakpoints.}

\author{W. Huber \email{huber@ebi.ac.uk}}

\keyword{manip}
