\name{normalizeByReference}
\alias{normalizeByReference}

\title{
  Probe-response normalization of hybridization intensities from a short-oligo microarray
}
\description{
  Adjust the hybridization intensities from a short-oligo microarray
  for probe-specific response effect by using a set of reference hybridizations.
}
\usage{
normalizeByReference(x, reference, pm, background, nrStrata=10,
  cutoffQuantile=0.05, plotFileNames, verbose=FALSE)
}
\arguments{
  \item{x}{eSet containing the data to be normalized.}
  \item{reference}{eSet with the same number of features as \code{x},
    containing the reference signal.}
  \item{pm}{indices specifying the perfect match features in
    \code{reference}.
    This can be either an integer vector with values between 1 and
    \code{nrow(exprs(reference))} or a logical vector.}
  \item{background}{indices specifying a set of background features in
    \code{x}.
    This can be either an integer vector with values between 1 and
    \code{nrow(exprs(x))} or a logical vector.}
  \item{nrStrata}{Integer (length 1), number of strata for estimating of
    additive background function.}
  \item{cutoffQuantile}{Numeric (length 1), the probes whose reference
    signal is below this quantile are thrown out.}
  \item{plotFileNames}{Character vector whose length is the same as the
    number of arrays in \code{x}. Optional, if missing, no plots are
    produced.}
  \item{verbose}{Logical of length 1, if \code{TRUE}, some messages about
    progress are printed.}
}
\details{...}

\value{
  An object of class \code{\link{eSet}}.
}

\author{W. Huber \email{huber@ebi.ac.uk}}

\references{...}

\examples{
  ##
}
\keyword{manip}
