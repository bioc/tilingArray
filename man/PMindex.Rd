\name{PMindex}
\alias{PMindex}
\alias{BGindex}

\title{Find the index of the exact match (PM) or background probes from a probeAnno environment}
\description{Find the index of the exact match (PM) or background probes from a probeAnno environment}
\usage{
PMindex(probeAnno)
BGindex(probeAnno)
}
\arguments{
\item{probeAnno}{environment with probe annotations. See
    package \code{davidTiling} for an example (?probeAnno).}
}
\details{
  These functions extract the exact match probes (PM) or background
  probes (from intergenic regions outside of known annotations) indices from \code{probeAnno}.  
  These indices can be used to select the relevant rows of intensity
  data from the \code{ExpressionSet} object for plotting and normalization.
}

\value{Numeric vector of indices.}

\author{Matt Ritchie <ritchie@ebi.ac.uk>}

\examples{
## library(davidTiling)
## data(davidTiling)
## data(probeAnno)
## pmind <- PMindex(probeAnno)
## mmind <- MMindex(probeAnno)
## bgind <- BGindex(probeAnno)
## boxplot(as.data.frame(log2(exprs(davidTiling))[pmind,]), outline=FALSE)           
}
\keyword{manip}
