\name{normalizeByReference}
\alias{normalizeByReference}

\title{
  Probe-specific normalization of hybridization intensities from an oligonucleotide microarray
}
\description{
  Adjust the hybridization intensities from an oligonucleotide microarray
  for probe-specific response effect by using one or several reference
  hybridizations.
  If \code{x} contains more than one array, \code{\link[vsn:vsn2]{vsnMatrix}} from
  the \code{vsn} package is called for between array normalization.
}
\usage{
normalizeByReference(x, reference, pm, background, refSig, nrStrata=10,
  cutoffQuantile=0.05, plotFileNames, verbose=FALSE)
}
\arguments{
  \item{x}{ExpressionSet containing the data to be normalized.}
  \item{reference}{ExpressionSet with the same number of features as \code{x},
    containing the reference signal, on the raw scale (non-logarithmic).
    This argument can be used to directly input the data from a set of
    replicate DNA hybridizations.
    Alternatively, the argument \code{refSig} can be specified.}
  \item{pm}{Indices specifying the perfect match features in
    \code{reference} (see Details).
    This can be either an integer vector with values between 1 and
    \code{nrow(exprs(reference))} or a logical vector.}
  \item{background}{Indices specifying a set of background features in
    \code{x} (see Details).
    This can be either an integer vector with values between 1 and
    \code{nrow(exprs(x))} or a logical vector.}
  \item{refSig}{A numeric vector of the same length as \code{pm} with
    estimates of probe response effects, on a logarithm-like scale.
    This argument can be specified alternatively to \code{reference}.}
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
  reference values.
  Typically, the reference values are obtained by hybridizing a DNA
  sample to the array, so that the abundance of target is
  the same for all reference probes, and their signal can be used to estimate
  the probe sequence effect. A reference probe is a probe that
  perfectly matches the target genome exactly once.
  Usually, not all probes on a chip are reference
  probes, hence the subset of those that are is specified by the
  argument \code{pm}.

  The background signal is estimated from the probes indicated by the
  argument \code{background}. They need to be a strict subset of the
  \code{reference} probes. I.e., they need to uniquely
  match the target organism's DNA, but are not expected to match any of its transcripts. A
  robust estimation method is used, so a small fraction of
  \code{background} probes that do hit transcripts is not harmful.

  A limitation of this normalization method is that it only makes sense
  for the data from reference probes, NA values are returned for all other probes.

  The functions \code{PMindex} and \code{BGindex} can be used to produce
  the \code{pm} and \code{background} arguments from a \code{probeAnno}
  environment such as provided in the davidTiling package.

  To summarize, a reference probe (indicated by argument \code{pm}) is a
  probe that perfectly matches the target genome exactly once, a
  background probe (indicated by argument \code{background})
  is a reference probe which we expect not to be transcribed. These
  should not be confused with what is called 'perfect match' and
  'mismatch' probes in Affymetrix annotation.
}

\seealso{\code{PMindex}, \code{BGindex}}

\value{
  A copy of \code{x} with the normalized intensities.
}

\author{W. Huber \email{huber@ebi.ac.uk}}

\references{The method implemented in this function is described in
  detail in Section 2.3 of the article
  Huber W, Toedling J, Steinmetz, L. Transcript mapping with
  high-density oligonucleotide tiling arrays. Bioinformatics 22, 1963-1970 (2006).}
\examples{
  ## see vignette assessNorm.Rnw in inst/scripts directory
}
\keyword{manip}
