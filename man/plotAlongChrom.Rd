\name{plotAlongChrom}
\alias{plotAlongChrom}
\title{Plot the results of segmentation along a chromosomal region}
\keyword{hplot}
\description{Plot the results of segmentation along a chromosomal region}

\usage{
plotAlongChrom(segObj, y, probeAnno, gff,
    nrBasesPerSeg,
    isDirectHybe=FALSE, 
    what = c("dots"), ## "heatmap"
    chr, coord, highlight, ylim, 
    colors, 
    doLegend=TRUE,
    featureColorScheme=1,
    featureExclude=c("chromosome","gene","nucleotide_match", "insertion", "intron"),
    featureNoLabel=c("CDS", "uORF"),
    main, ...)
}

\arguments{
  \item{segObj}{environment containing three types of data
    \itemize{
      \item microarray intensities in along-chromosome order,
      \item a segmentation (output of findSegments),
      \item a data.frame named \code{segScore} with segment scores; can
        be missing iff \code{nrBasesPerSeg} is present,
      \item a numeric scalar names \code{theThreshold}, which is used to 
        draw a horizontal "threshold" line in the plot.
    }
    See also Paragraph \emph{Intensities} in the \emph{Details} section.
    This is currently simply an environment that is expected to contain
    certain objects with certain names. In future version, this is
    may be replaced by a proper S4 object.}
  \item{y}{a numeric vector or matrix containing the signal to be
    plotted. 
    See also Paragraph \emph{Intensities} in the \emph{Details} section.}
  \item{chr}{integer of length 1 with the number of the chromosome for
    which to make the plot.}
  \item{coord}{integer of length 2 with start and end coordinates of the
    plot.}
  \item{highlight}{optional, list with two elements: a single numeric value
    \code{coord} and a character \code{strand}. If present, this
    position is marked by a vertical red bar on the coordinate axis.}
  \item{ylim}{numeric vector of length two with y-axis limits.}
  \item{nrBasesPerSeg}{numeric scalar. If \code{segObj} does not 
     contain a \code{segScore} data.frame, then this parameter allows to 
     choose the parameter which controls the number of segments in the 
     segmentation (called \emph{S} in the paper). In particular, this 
     parameter specifies the average length of segments, measured in bases, 
     thus one value can be used for all chromosomes.}
  \item{probeAnno}{environment with probe annotations, see vignette.}
  \item{gff}{data frame with genome annotation from the GFF file.}
  \item{colors}{named character vector, optional. If missing, 
    the following default is used:
    \code{c("+"="#00441b", "-"="#081d58", "duplicated"="grey", "cp"="#101010", 
   "highlight"="red", "threshold"="grey")}, 
   where the first three elements refer to colors of data points and the 
   last three to those of lines in the plot.}
  \item{featColScheme}{numeric scalar, used to select a color scheme for the 
   boxes representing genomic features such as coding sequences, ncRNAs etc. 
   Currently the values 1 and 2 are supported.}
  \item{isDirectHybe}{logical scalar: if TRUE, the mapping of probes to genomic 
   strands is reversed with respect to the default. This is appropriate for data 
   from a direct RNA hybridization that used no reverse transcription.} 
  \item{scoreShow}{character: name of a column in \code{segScore} that
    is to be used for the coloring of the segments.}
  %\item{haveNames}{logical; should their names be added to the feature symbols?}
  \item{haveLegend}{logical: should the plot contain a legend?} 
  %\item{main}{character: plot title.}
  %\item{pointSize}{unit object: point size used for the probe intensities scatterplot.}
  \item{...}{further arguments passed on to \code{plot.default}}
}

\details{
  \emph{Intensities}: There are two alternative, mutually exclusive
  ways of providing the intensities to be plotted to this function.
  \enumerate{
    \item Via the parameters \code{y} and \code{probeAnno}. Then,
    \code{y} is a vector or matrix of intensities.
    If it is a matrix, its rows correspond to probes on the array, and
    its columns to different conditions, time points, etc. 
    \code{probeAnno} is an
    environment that contains integer vectors of the same length as
    \code{y} (if it is a vector)

    with start positions along
    the chromosomes and indices in \code{y} of the probes. For example,
    the start positions and indices of probes for the + strand of
    chromosome 1 would be described by environment elements
    \code{"1.+.start"} and \code{"1.+.index"}, respectively.
    \item Via the parameter \code{segRes}. In that case,
    the intended workflow for using this function is as follows:
    First, use the script \code{segment.R} to generate a segmentation.
    This can be run in parallel on several processors, separately for each
    chromosome and strand. The results of this are stored in files of the
    name \code{1.+.rda}, \code{1.-.rda}, \code{2.+.rda}, and so forth,
    typically within a dedicated directory.
    Then, the script \code{viewsegmentation.R} can be used to view the
    segmentation results. It collects the \code{R} objects in these
    \code{.rda} files into the environment \code{segRes}, with which it calls 
    this function.

    Furthermore, the script \code{scoreSegments.R} can be used to
    calculate scores for each segment, which can be passed to this
    function using the \code{segRes} argument.
  }
} 
  

\author{Wolfgang Huber <huber@ebi.ac.uk>}
\examples{
  ## 
 }

