\name{plotAlongChrom2}
\alias{plotAlongChrom2}
\title{Plot the results of segmentation along a chromosomal region}
\keyword{hplot}
\description{Plot the results of segmentation along a chromosomal region}

\usage{plotAlongChrom2(chr, coord, highlight, segRes, segScore,
   scoreShow = "pt", nrBasesPerSeg, gff)}

\arguments{
  \item{chr}{integer of length 1 with the number of the chromosome for
    which to make the plot.}
  \item{coord}{integer of length 2 with start and end coordinates of the
    plot.}
  \item{highlight}{optional, list with two elements: a single numeric value
    \code{coord} and a character \code{strand}. If present, this
    position is marked by a vertical red bar on the coordinate axis.}
  \item{segRes}{environment containing microarray intensities in
    along-chromosome order and segmentation results, see details.}
  \item{segScore}{dataframe with segmentation scores, see details. Can
    be missing iff \code{nrBasesPerSeg} is present.}
  \item{scoreShow}{character, name of a column in \code{segScore} that
    is to be used for the coloring of the segments.}
  \item{nrBasesPerSeg}{numeric of length 1 with the average number of
    nucleotide bases per segment. Can
    be missing iff \code{segScore} is present.}
  \item{gff}{data frame with genome annotation from the GFF file.}
}

\details{The intended workflow for using this function is as follows:
  First, use the script \code{segment.R} to generate a segmentation.
  This can be run in parallel on several processors, separately for each
  chromosome and strand. The results of this are stored in files of the
  name \code{1.+.rda}, \code{1.-.rda}, \code{2.+.rda}, and so forth,
  typically within a dedicated directory.
  Then, the script \code{viewsegmentation.R} can be used to view the
  segmentation results. It collects the \code{R} objects in these
  \code{.rda} files into the environment \code{segRes}. It then calls 
  this function.
  
  Furthermore, the script \code{scoreSegments.R} can be used to
  calculate scores for each segment, which can be passed tpo this
  function using the \code{segRes} argument.}

\author{Wolfgang Huber <huber@ebi.ac.uk>}
\examples{
  ## 
 }

