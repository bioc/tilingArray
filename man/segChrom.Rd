\name{segChrom}
\alias{segChrom}

\title{Fit a piecewise constant curve to along chromosome data (wrapper function)}
\description{
  Wrapper around the \code{segment} function for each strand of one or
  more chromosomes specified by the user. It does some typical
  preprocessing and I/O.
}
  
\usage{
segChrom(y, probeAnno, chr=1:17, strands=c("+", "-"), 
  nrBasesPerSegment = 1500, maxk = 3000, step = 7, confint = FALSE,
  confintLevel = 0.95, useLocks=TRUE, verbose=TRUE, savedir)
}

\arguments{  
  \item{y}{ExpressionSet or matrix containing the data to be segmented.}
  \item{probeAnno}{environment with probe annotations. See
    package \code{davidTiling} for an example (?probeAnno).}
  \item{chr}{integer scalar or vector specifying which chromosome(s) to segment.}
  \item{strands}{character scalar or vector specifying which strands to segment.}
  \item{nrBasesPerSegment}{integer (length 1): the parameter
    \code{maxseg} of the \code{\link{segment}} function is calculated as the length of the chromosome divided
    by \code{nrBasesPerSegment}. Thus, it determines the average segment
    length in the finest segmentation.}
  \item{maxk}{passed on to  the function \code{\link{segment}}.}
  \item{step}{integer scalar, indicating the minimum distance between
    consecutive probes. In cases when probes are offset by less than
    \code{step} bases, the probes are sampled to achieve the desired spacing.}
  \item{confint}{logical scalar.  If \code{TRUE}, confidence intervals
    for each change-point are calculated.}
  \item{confintLevel}{numeric scalar between 0 and 1 indicating the
    probability level for the confidence intervals that are calculated
    for each change-point.}
  \item{useLocks}{logical scalar. Should a file locking mechanism be used
    that allows for a simple-minded parallelization of this function.}
  \item{verbose}{logical scalar. Should we be chatty about our progress?}
  \item{savedir}{character scalar. If specified, resulting \code{segmentation}
    objects are saved (with \code{\link{save}}) to this directory.}

}

\details{
  This function is a wrapper for the \code{\link{segment}} function.
  Refer to its help page for further details.
}

\value{
  An environment containing S4 objects of class \code{"segmentation"}
  called "1.+", "1.-", etc. (depending on the values in \code{chr} and
  \code{strands}),
  where "+" and "-" indicate the strand and the preceding number refers
  to the chromosome.
  If \code{savedir} is specified, there is also the side-effect that
  a series of files "1.+.rda", "1.-.rda", etc. is saved in that directory.}

\author{Matt Ritchie <ritchie@ebi.ac.uk> and Wolfgang Huber <huber@ebi.ac.uk>}

\examples{
\dontrun{
  library("davidTiling")
  data("davidTiling")
  data("probeAnno")
  isDNA = seq(1:3)
  yn = normalizeByReference(davidTiling[,-isDNA],davidTiling[,isDNA], probeAnno=probeAnno)
  seg = segChrom(yn, probeAnno) ## this will take a while to run!
}
}
\keyword{manip}
