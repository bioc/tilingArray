\name{scoreSegments}
\alias{scoreSegments}
\title{Score segments}
\description{Score the segments found by a previous call to findSegments by
  comparing to genome annotation}
\usage{
scoreSegments(s, gff, 
  nrBasePerSeg = 1500, 
  probeLength  = 25,
  params = c(overlapFraction = 0.5, oppositeWindow = 100, flankProbes=10),
  verbose = TRUE)
}
\arguments{
  \item{s}{environment. See details.}
  \item{gff}{GFF dataframe.}
  \item{nrBasePerSeg}{Numeric of length 1. This parameter determines the
    number of segments.}
  \item{probeLength}{Numeric of length 1.}
  \item{params}{vector of additional parameters, see details.}
  \item{verbose}{Logical.}
}

\details{This function scores segments.
It is typically called after a \emph{segmentation}.
For an example segmentation script, see the script \code{segment.R}
in the \code{scripts} directory of this package. For an example
scoring script, which loads the data and then calls this function, see
the script \code{scoreSegments.R}.

More details can be found in the separate file mat_meth_wh_DATE.doc.


% For each segment, we calculate and record its:
% \describe{
% \item{chr, strand}{chromosome and strand}
% \item{start, end, length}{start position, end position, length (in
% bases)}
% \item{frac.dup}{fraction (0...1) of probes in this segment that have
% also hits otherwhere} 
% \item{level}{mean signal level}
% \item{geneInSegment}{list of gene identifiers (these can be zero, one,
%   or several identifiers). A gene is included if it is fully contained
%   within the segment, i.e. its start coordinates are >= the segment's
%   start and its end coordinates <= the segment's end.}
% \item{overlappingFeature}{list of feature identifiers (these include
%   genes, CDSs (=exons), ncRNAs ..., everything which has a line the
%   GFF file). A feature is included if the overlap between it and the
%   segment is more than a fraction \code{minOverlapFractionSame} of
%   the length of the segment or of the feature, whichever is smaller. The
%   overlappingFeature list by definition contains the elements of the
%   geneInSegment list, but can be larger.}
% \item{oppositeFeature}{list of feature identifiers. A feature is
%   included if the overlap between it and the segment is >=
%   minOverlapOppo (see below) bases.}
% \item{oppositeExpression}{a number. The signal on the opposite strand is
%   filtered with a moving average smoother of width oppositeWindow (see
%   below) bases. oppositeExpression is the minimum of the result. It is
%   later used to eliminate potential reverse transcription artifacts from
%   the unannotated, potential antisense segments. If it is
%   sufficiently small, we can assume that there is no transcription at
%   least on parts of the strand opposite the segment, hence it cannot be
%   a reverse transcription artifact.}
% \item{isIsolatedSame, isIsolatedOppo}{logical. TRUE if the distance
%   between the segment and and any annotated feature on the same or
%   opposite strand, respectively, is >=
%   \code{params["minIsolatedDistance"]} bases.}
% \item{utr3, utr5}{positive integer numbers. These are calculated only
%   for segments which have exactly one gene in the geneInSegment
%   list. They are calculated as the difference between start points of
%   segment and gene, and between end points of segment and gene,
%   respectively.}
% \item{distLeft, distRight}{positive integer numbers. distLeft is the
%   distance between the start of the segment and the closest end of any
%   annotated feature, and distRight is the distance between the end of
%   the segment and the closest start any annotated feature.}
% \item{zLeft, zRight}{z-scores for the left (right) flank of the
%   segment. They are calculated as the difference between mean and of the
%   segment and mean of the signal from a region of length utrScoreWidth
%   (see below) immediately to the left (right), divided by the standard
%   deviation of the region. Note that the standard deviation of the
%   signal within the segment is not considered here.}
% }

% The meaning of the parameters in the parameter vector \code{params} is
% as follows:
% \describe{
%   \item{minOverlapFractionSame}{see the definition of
%     \code{overlappingFeature} above}
%   \item{minOverlapOppo}{see the definition of
%     \code{oppositeFeature} above}
%   \item{oppositeWindow}{see the definition of \code{oppositeExpression}
%     above}
%   \item{utrScoreWidth}{see the definition of \code{zLeft}, \code{zRight}
%     above}
% }
}


\value{
  A dataframe with columns as described in the details section.
}

\author{W. Huber <huber@ebi.ac.uk>}

\examples{
}
\keyword{manip}
