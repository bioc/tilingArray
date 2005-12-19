\name{segmentation-class}
\docType{class}
\alias{segmentation-class}
\alias{confint,segmentation-method}
\alias{plot,segmentation-method}
\alias{show,segmentation-method}

\title{Class "segmentation" representing a segmentation result}
\description{This class represents the result of a segmentation, usually
a call of the function \code{findSegments} or \code{segment}.}
\section{Objects from the Class}{
  Objects can be created by calls of the function
  \code{findSegments} or \code{segment}.
}
\section{Slots}{
	 \describe{
    \item{\code{y}:}{The numeric data that has been segmented.}
    \item{\code{breakpoints}:}{Numeric vector of change-points,
      i.e. positions, at which a new segment starts in the
      numeric data \code{y}}
    \item{\code{negloglik}:}{negative log-likelihood of the segment
      model under the data \code{y}}
    \item{\code{hasConfint}:}{logical; have confidence intervals been
      computed for the change-points.}
  }
}
\section{Methods}{
  \describe{
    \item{confint}{\code{signature(object = "segmentation")}: compute
      confidence intervals for the change-points of the segmentation.}
    \item{plot}{\code{signature(x = "segmentation")}: plot the result of
    segmentation, i.e. the numeric data along with segment borders and
    possibly confidence intervals for these as well.}
    \item{show}{\code{signature(object = "segmentation")}: show summary
      of the segmentation object}
  }
}
%\references{ ~put references to the literature/web site here ~ }
\author{Wolfgang Huber \email{huber@ebi.ac.uk}, Joern Toedling
  \email{toedling@ebi.ac.uk}}
%\note{ ~~further notes~~ }
\seealso{\code{\link{findSegments}}}
%\examples{}
\keyword{classes}
