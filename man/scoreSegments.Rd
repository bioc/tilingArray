\name{scoreSegments}
\alias{scoreSegments}
\title{Score segments}
\description{Score the segments found by a previous call to findSegments by
  comparing to genome annotation}
\usage{
scoreSegments(s, gff, 
  nrBasePerSeg = 1500, 
  probeLength  = 25,
  knownFeatures = c("CDS", "gene", "ncRNA", "nc_primary_transcript",
        "rRNA", "snRNA", "snoRNA", "tRNA",
        "transposable_element", "transposable_element_gene"),
  params =c(minOverlapFractionSame=0.8, minOverlapBasesOppo=40), 
  verbose = TRUE)
}
\arguments{
  \item{s}{environment. See details.}
  \item{gff}{GFF dataframe.}
  \item{nrBasePerSeg}{Numeric of length 1. This parameter determines the
    number of segments.}
  \item{probeLength}{Numeric of length 1.}
  \item{knownFeatures}{Character vector. Names of those features in
    \code{gff} which should be considered \emph{known features}.}
  \item{params}{vector of additional parameters, describing limits for
    probes considered for plotting}
  \item{verbose}{Logical.}
}

\details{
This function is typically called after a \emph{segmentation}.
For an example segmentation script, see the script \code{segment.R}
in the \code{scripts} directory of this package. For an example
scoring script, which loads the data and then calls this function, see
the script \code{scoreSegments.R}.
  
This function scores segments.
For each segment, we calculate and/or record its:
\describe{
\item{chr, strand}{chromosome and strand}
\item{start, end, length}{start position, end position, length (in
bases)}
\item{frac.dup}{fraction (0...1) of probes in this segment that have
also hits otherwhere} 
\item{level}{mean level minus that of the background}
\item{featureInSegment}
{the name of one or several features that are contained to more
than \code{minOverlapFractionSame} in the segment. Can be empty.}
\item{segmentInFeature}
{the name of the (single) feature that the segment is contained in
(to 100 percent). Can be empty. This is used in the UTR mapping.}
\item{oppositeFeature}
{the name of one several features that have >=\code{minOverlapBasesOppo} bases
overlap with the segment.}
\item{utr5}{candidate 5' UTR length}
\item{utr3}{candidate 5' UTR length}
\item{and other}{...to be filled in...}
}
A probe is considered \emph{part of} a feature if its middle
nucleotide (\code{(probeLength-1)/2}) falls into it. That is an
arbitrary rule, which is intended to err on the side of assigning
probes to features rather than not.
}


\value{
  A dataframe with columns as described in the details section.
}

\author{W. Huber <huber@ebi.ac.uk>}

\examples{
}
\keyword{manip}
