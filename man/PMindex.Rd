\name{PMindex, MMindex, BGindex}
\alias{PMindex}
\alias{MMindex}
\alias{BGindex}

\title{Find the index of the PM, MM or background probes from a probeAnno environment}
\description{Find the index of the PM, MM or background probes from a probeAnno environment}
\usage{
PMindex(probeAnno)
MMindex(probeAnno)
BGindex(probeAnno)
}
\arguments{
\item{probeAnno}{environment with probe annotations. See
    package \code{davidTiling} for an example (?probeAnno).}
}
\details{
  These functions extract the perfect match (PM), mismatch (MM) or background probe (from intergenic regions outside of known annotations) indicies from the probeAnno.  
These indices can be used to select the relevant rows of intensity data from the eSet object (read in from the cel files) for plotting and normalization.
}

\value{Numeric vector returning of indicies.
}

\author{Matt Ritchie <ritchie@ebi.ac.uk>}

\examples{## library(davidTiling)
## data(davidTiling)
## data(probeAnno)
## pmind <- PMindex(probeAnno)
## mmind <- MMindex(probeAnno)
## bgind <- BGindex(probeAnno)
## boxplot(as.data.frame(log2(exprs(davidTiling))[pmind,]), outline=FALSE)           
}
\keyword{manip}
