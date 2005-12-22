\name{segmentation-class}
\docType{class}
\alias{segmentation}
\alias{class:segmentation}


\concept{confint}
\concept{plot}
\concept{show}

\title{
  The class segmentation represents a segmentation result.
}
\description{
  This class represents the result of a segmentation, usually
  a call to the function \code{segment}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the function
  \code{segment} or by calls of the form \code{new("segmentation", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{y}:}{A matrix with the data. Rows correspond to the
      \code{x}-variable, columns to replicate measurements at the same
      value of \code{x}. Breakpoints are fitted along the \code{x}-axis.
      For example \code{x} can be genomic coordinates or time.} 
    \item{\code{breakpoints}:}{List of segmentations. The element
      \code{breakpoints[[j]]} corresponds to a segmentation into
      \code{j} segments, i.e. with \code{j-1} breakpoints. It is a
      matrix with \code{(j-1)} rows and 1 or 3 columns. It always
      contains a columns named \code{estimate} with the point estimates.
      Optionally, it may contain columns \code{lower} and \code{upper}
      with the confidence intervals. The point estimates are the row
      indices in \code{y} where new segments start, for example:
      let \code{z=breakpoints[[j]]}, then the first segment
      is from row \code{1} to \code{z[1, "estimate"]},
      the second from row \code{z[1, "estimate"]} to
      \code{z[2, "estimate"]-1}, and so on.}
    \item{\code{negloglik}:}{numeric vector of the same length as
      \code{breakpoints}. The
      negative log-likelihood of the piecewise
      constant models under the data \code{y}.}
    \item{\code{hasConfint}:}{logical vector of the same length as
      \code{breakpoints}. TRUE if the confidence interval estimates
      are present, i.e. if the matrix \code{breakpoints[[j]]} has
      columns \code{lower} and \code{upper}.}
  }
}

\section{Methods}{
  \describe{
    \item{confint}{The method \code{confint(object, parm, level=0.95,
	het.reg=FALSE, het.err=FALSE, ...)} computes confidence
      intervals for the change point estimates of the
      segmentation. Typically, these were obtained from a previous call
      to the function \code{link{segment}} that created the object.
      This is just a wrapper for the function
      \code{\link[strucchange]{confint.breakpointsfull}}
      from the \code{strucchange} package, which does all the hard
      computations.
      Parameters: \code{object} an object of class \code{segmentation},
      \code{parm} an integer vector, it determines for which of the segmentation fits
      confidence intervals are computed. See also \code{\link{segment}}.
      The other parameters are directly passed on to
      \code{\link[strucchange]{confint.breakpointsfull}}.
    }
    \item{plot}{The method \code{plot(x, y, xlim, xlab="x", ylab="y", bpcol, pch=16, ...)}
      provides a simple visualization of the result of a
      segmentation. Parameters: \code{x} an object of class \code{segmentation},
      \code{y} an integer between \code{1} and
      \code{length(x@breakpoints)}, selecting which of the fits
      contained in \code{x} to plot. The plot shows the numeric data along with segment borders and
      if available confidence intervals.}
    \item{show}{summary.}
  }
}

%\references{ ~put references to the literature/web site here ~ }
\author{Wolfgang Huber \email{huber@ebi.ac.uk}, Joern Toedling
  \email{toedling@ebi.ac.uk}}

\seealso{\code{\link{segment}}}
\keyword{classes}

\example{
## generate random data with 5 segments:
x = unlist(lapply(c(0,3,0.5,1.5,5), function(m) rnorm(10, mean=m)))

seg = segment(x, maxseg=10, maxk=15)
seg = confint(seg, nSegments=c(3,4,5,6))

# see that the step between segments 3 and 4 is less certain than
#  the other ones:
if(interactive())
  plot(seg)

show(segres)
}

