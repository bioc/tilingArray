\name{scoreSegments}
\alias{scoreSegments}
\title{Score segments}
\description{Score the segments found by a previous call to findSegments by
  comparing to genome annotation}
\usage{
scoreSegments(x, gff, 
  nrBasePerSeg = 1500, 
  probeLength  = 25,
  knownFeatures = c("CDS", "gene", "ncRNA", "nc_primary_transcript",
        "rRNA", "snRNA", "snoRNA", "tRNA",
        "transposable_element", "transposable_element_gene"),
  verbose = TRUE)
}
\arguments{
  \item{x}{environment. See details.}
  \item{gff}{GFF dataframe.}
  \item{nrBasePerSeg}{Numeric of length 1. This parameter determines the
    number of segments.}
  \item{probeLength}{Numeric of length 1.}
  \item{knownFeatures}{Character vector. Names of those features in
    \code{gff} which should be considered \emph{known features}.}
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
    \item{p-value}{one-sample t-test against hypothesis level=0}
    \item{same.feature}{the name(s) of annotated features on the same
      strand with which it overlaps}
    \item{same.overlap}{the fraction (0...1) of probes that are part of
      an annotated feature (see below for definition of \emph{part of})}
    \item{same.dist.start2feat}{Here, we distinguish two cases:
      \enumerate{
         \item If the leftmost probe of a segment is part of an annotated feature,
          then dist.start2feat is the difference between that probe's middle
          nucleotide and the annotated start point of the feature.
         \item If it is not part of an annotated feature, then dist.start2feat
          is the difference between the middle nucleotide and the end point of
          the next feature to the left.}}
    \item{same.dist.end2feat}{Here, we distinguish two cases:
      \enumerate{
         \item If the rightmost probe of a segment is part of an annotated feature,
          then dist.endfeat is the difference between that probe's middle
          nucleotide and the annotated end point of the feature.
         \item If it is not part of an annotated feature, then dist.end2feat
          is the difference between the start point of the next feature to the
          right and the middle nucleotide.}}}
    In addition, quantities oppo.XXX are calculated for the opposite
    strand analogous to the same.XXX.

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
