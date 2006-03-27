\name{segmentation}
\docType{class}
\alias{segmentation-class}

\alias{confint}
\alias{logLik}
\alias{confint,segmentation-method}
\alias{logLik,segmentation-method}
\alias{plot,segmentation-method}
\alias{show,segmentation-method}

\concept{plot}
\concept{show}
\concept{confint}
\concept{logLik}

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
    \item{\code{y}:}{A matrix with the data (the dependent variable(s)), see \code{\link{segment}}.} 
    \item{\code{x}:}{A numeric vector with the regressor variable.
      The length of this vector must be either the same as
      \code{nrow(y)}, or 0. The latter case is equivalent to
      \code{x=1:nrow(y)}.}
    \item{\code{flag}:}{An integer vector, whose length must be either the same as
      \code{nrow(y)}, or 0. This can be used to \emph{flag} certain
      probes for special treatment, for example by
      \code{\link{plotAlongChrom}}.}
    \item{\code{breakpoints}:}{List of segmentations. The element
      \code{breakpoints[[j]]} corresponds to a segmentation fit of
      \code{j} segments, i.e. with \code{j-1} breakpoints. It is a
      matrix with \code{(j-1)} rows and 1 or 3 columns. It always
      contains a column named \code{estimate} with the point estimates.
      Optionally, it may contain columns \code{lower} and \code{upper}
      with the confidence intervals. The point estimates are the row
      indices in \code{y} where new segments start, for example:
      let \code{z=breakpoints[[j]]}, then the first segment
      is from row \code{1} to \code{z[1, "estimate"]-1},
      the second from row \code{z[1, "estimate"]} to
      \code{z[2, "estimate"]-1}, and so on.}
    \item{\code{negloglik}:}{Numeric vector of the same length as
      \code{breakpoints}. The
      negative log-likelihood of the piecewise
      constant models under the data \code{y}.}
    \item{\code{hasConfint}:}{Logical vector of the same length as
      \code{breakpoints}. TRUE if the confidence interval estimates
      are present, i.e. if the matrix \code{breakpoints[[j]]} has
      columns \code{lower} and \code{upper}.}
    \item{\code{nrSegments}:}{A scalar integer, value must be either
      \code{NA} or between \code{1} and \code{length(breakpoints)}.
      Can be used to select one of the fits in \code{breakpoints} for
      special treatment, for example by
      \code{\link{plotAlongChrom}}.}
  }
}

\section{Methods}{
  \describe{
    \item{confint}{The method \code{confint(object, parm, level=0.95,
	het.reg=FALSE, het.err=FALSE, ...)} computes confidence
      intervals for the change point estimates of the
      segmentation. Typically, these were obtained from a previous call
      to the function \code{\link{segment}} that created the object.
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
    \item{logLik}{The method \code{logLik(object, penalty="none", ...)}
      returns the log-likelihoods of fitted models. Valid values for the argument
      \code{penalty} are \code{none}, \code{AIC} and \code{BIC}.}
    \item{plot}{The method \code{plot(x, y, xlim, xlab="x", ylab="y",
	bpcol="black", bplty=1, pch=16, ...)}
      provides a simple visualization of the result of a
      segmentation. Parameters: \code{x} an object of class \code{segmentation},
      \code{y} an integer between \code{1} and
      \code{length(x@breakpoints)}, selecting which of the fits
      contained in \code{x} to plot, \code{bpcol} and \code{bplty} color
      and line type of breakpoints. The plot shows the numeric data
      along with breakpoints and if available their confidence intervals.}
    \item{show}{summary.}
  }
}

%\references{ ~put references to the literature/web site here ~ }
\author{Wolfgang Huber \email{huber@ebi.ac.uk}}

\seealso{\code{\link{segment}}}
\keyword{classes}

\examples{
## generate random data with 5 segments:
y = unlist(lapply(c(0,3,0.5,1.5,5), function(m) rnorm(10, mean=m)))

seg = segment(y, maxseg=10, maxk=15)
seg = confint(seg, parm=c(3,4,5))

if(interactive())
  plot(seg, 5)

show(seg)
}

