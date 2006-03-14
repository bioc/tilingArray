\name{plotAlongChrom}
\alias{plotAlongChrom}
\title{Plot signals and segmentation along a chromosomal region}
\keyword{hplot}
\description{Plot signals and segmentation along a chromosomal region}

\usage{
plotAlongChrom(segObj, y, probeAnno, gff,
    isDirectHybe=FALSE, 
    what = c("dots"), ## "heatmap"
    chr, coord, highlight, ylim, 
    colors, 
    doLegend=TRUE,
    featureColorScheme=1,
    featureExclude=c("chromosome","gene","nucleotide_match", "insertion", "intron"),
    featureNoLabel=c("uORF"),
    pointSize=unit(0.6, "mm"),
    main, ...)
}

\arguments{
  \item{segObj}{Either an environment or an object of S4 class
    \code{segmentation}.
    See \emph{Details}.}
  \item{y}{a numeric vector or matrix containing the signal to be
    plotted. See \emph{Details}.}
  \item{probeAnno}{environment with probe annotations. See
    \emph{Details}, and package
    \code{davidTiling} for an example.}
  \item{gff}{data frame with genome annotation from the GFF file.}
  \item{isDirectHybe}{logical scalar: if TRUE, the mapping of probes to genomic 
   strands is reversed with respect to the default. This is appropriate for data 
   from a direct RNA hybridization that used no reverse transcription.} 
  \item{what}{character scalar, choice of signal visualization, can be
    either \code{dots} or \code{heatmap}}  
  \item{chr}{integer of length 1 with the number of the chromosome for
    which to produce the plot.}
  \item{coord}{integer of length 2 with start and end coordinates (in
    bases) for which to produce the plot.}
  \item{highlight}{optional, list with two elements: a single numeric value
    \code{coord} and a character \code{strand}. If present, this
    position is marked by a vertical red bar on the coordinate axis.}
  \item{ylim}{numeric vector of length two with y-axis limits.}
  \item{colors}{named character vector, optional. If missing, 
    a default color scheme is used:
    \code{c("+"="#00441b", "-"="#081d58", "duplicated"="grey", "cp"="#101010", 
   "highlight"="red", "threshold"="grey")}, 
   where the first three elements refer to colors of data points and the 
   last three to those of lines in the plot.}
  \item{doLegend}{logical: should the plot contain a legend?} 
  \item{featureColorScheme}{numeric scalar, used to select a color scheme for the 
   boxes representing genomic features such as coding sequences, ncRNAs etc. 
   Currently the values 1 and 2 are supported.}
  \item{featureExclude}{character vector, names of feature types (in
    gff) that should not be plotted.
    Additional possible candidates include: "ARS", "repeat\_region",
    "repeat\_family", "nc\_primary\_transcript".
  }
  \item{featureNoLabel}{character vector, names of feature types (in
    gff) that should not be labeled with their names (if they are plotted).}
  \item{pointSize}{unit object: point size used for the probe intensities scatterplot.}
  \item{main}{character: plot title.}
  \item{...}{further arguments that can be passed on to the
    the functions that implement the \code{what} option (see above),
    \code{plotSegmentationDots} and \code{plotSegmentationHeatmap}.}
}

\details{
  \emph{Intensities}: There are two alternative, mutually exclusive
  ways of providing the intensities that are to be plotted to this function.
  \enumerate{
    \item Via the parameters \code{y} and \code{probeAnno}. In this case,
    \code{y} is a matrix of intensities, whose rows correspond to probes
    on the array, and its columns to different conditions, time points, etc.
    It is also acceptable that \code{y} is provided as a vector, in
    which case it is converted to an \code{nrow(y) x 1} matrix.
    \code{probeAnno} is an
    environment whose elements correspond to target sequences (e.g.
    chromosome strands) and that contain integer vectors of length
    \code{nrow(y)} with information about the probes: start and end positions of
    their alignment to the target sequence, their row indices in
    \code{y}, the type of alignment (is it perfect? is is unique?).
    For example,
    the start positions and indices of probes for the + strand of
    chromosome 1 would be described by environment elements
    \code{"1.+.start"} and \code{"1.+.index"}.
    
    \item Via the parameter \code{segObj}.
  }

  \emph{segObj}: This can be either an object of S4 class
  \code{segmentation} or an environment that by convention contains a
  certain set of objects.
  Future work on this package
  will focus on the S4 class \code{segmentation}. The \code{environment}
  option is provided for backward compatibility.
  
  \emph{Explanation of the environment}: the intended workflow is as follows:
  Use the script \code{segment.R} (in the \code{inst/scripts} directory of this
  package) to generate segmentations.
  This can be run in parallel on several processors, separately for each
  chromosome and strand. The results of this are stored in files of the
  name \code{1.+.rda}, \code{1.-.rda}, \code{2.+.rda}, and so forth,
  typically within a dedicated directory.
  Then use the script \code{readSegments.R} to collect the
  \code{R} objects in these \code{.rda} files into the environment.
  It contains three types of data:
  \itemize{
    \item microarray intensities in along-chromosome order.
    \item the segmentation objects (output of findSegments).
    \item a dataframe named \code{segScore} with segment scores; it can
    be missing iff \code{nrBasesPerSeg} is present,
    \item a numeric scalar names \code{theThreshold}, which is used to 
    draw a horizontal "threshold" line in the plot.
  }
} 
  

\author{Wolfgang Huber <huber@ebi.ac.uk>}
\examples{
  ## 1. see viewSegmentation.R script in the inst/scripts directory
  ## 2. (newer): segmentation.Rnw
 }

