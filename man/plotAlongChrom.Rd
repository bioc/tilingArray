\name{plotAlongChrom}
\alias{plotAlongChrom}
\title{Plot microarray data and annotation features along a chromosomal region}
\keyword{hplot}
\description{Plot microarray data and annotation features along a chromosomal region}

\usage{plotAlongChrom(y, chr, from, to, extend=0, gff, chrSeqname, probeAnno)}

\arguments{
  \item{y}{vector of real numbers, one for each probe on the array.}
  \item{chr}{character or integer of length 1: which chromosome?}
  \item{from}{integer of length 1: start coordinate}
  \item{to}{integer of length 1: end coordinate}
  \item{extend}{integer of length 1: extend start and end by how much}
  \item{gff}{data frame with genome annotation from the GFF file.}
  \item{chrSeqname}{character vector whose i-th element is the name of
    the i-th chromosome, as used in the \code{seqname} column of the
    \code{gff} table.}
  \item{probeAnno}{a probe annotation environment.}
}

\details{No details.}

\author{Wolfgang Huber <huber@ebi.ac.uk>}
\examples{
  ## TO DO
 }

