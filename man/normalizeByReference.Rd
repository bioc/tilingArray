\name{normalizeByReference}
\alias{normalizeByReference}

\title{
  Probe-specific normalization of hybridization intensities from a short-oligo microarray
}
\description{
  Adjust the hybridization intensities from a short-oligo microarray
  for probe-specific response effect by using one or several reference hybridizations.
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
  \item{nrStrata}{Integer (length 1), number of strata for the estimation of
    the background function.}
  \item{cutoffQuantile}{Numeric (length 1), the probes whose reference
    signal is below this quantile are thrown out.}
  \item{plotFileNames}{Character vector whose length is the same as the
    number of arrays in \code{x}. Optional, if missing, no plots are
    produced.}
  \item{verbose}{Logical of length 1, if \code{TRUE}, some messages about
    progress are printed.}
}
\details{The intensities in \code{x} are adjusted according to the
  reference values. For example, the reference values could be obtained
  by hybridizing a DNA sample to the array. The abundance of target is
  thus the same for all PM probes, and the intensities that are actually
  observed can be used for estimating the adjustment parameters. Since
  this only makes sense for the PM probes, NA values be returned for all
  other probes. The set of indices for background probes,
  \code{background}, is used for the
  estimation of the background correction.
  }

\value{
  An object of class \code{eSet}.
}

\author{W. Huber \email{huber@ebi.ac.uk}}

\references{Transcript mapping with 
  oligonucleotide high-density tiling arrays,
  W. Huber, J. Toedling, L. Steinmetz, preprint (2006).}
\examples{
  ## see vignette assessNorm.Rnw in inst/scripts directory
}
\keyword{manip}
