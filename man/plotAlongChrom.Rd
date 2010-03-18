\name{plotAlongChrom}
\alias{plotAlongChrom}
\title{Plot signals and segmentation for a region of a chromosome}
\keyword{hplot}
\description{Plot signals and segmentation for a region of a chromosome}

\usage{
plotAlongChrom(segObj, y, probeAnno, gff,
    isDirectHybe=FALSE, 
    what = c("dots"), ## "heatmap"
    chr, coord, highlight,
    colors, doLegend=FALSE,
    featureExclude=c("chromosome", "nucleotide_match", "insertion"),
    featureColorScheme=1, extras, 
    rowNamesHeatmap, rowNamesExtras, ylab, ylabExtras, main, 
    colHeatmap=colorRamp(brewer.pal(9, "YlGnBu")), 
    colExtras=colorRamp(brewer.pal(9, "Reds")), 
    sepPlots=FALSE, reOrder=TRUE, ...)
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
  \item{what}{character scalar indicating which signal visualization to plot.  Can be
    either \code{dots} to plot each probe intensity with a point, or \code{heatmap} to produce a colorscale representation of the intesities.}  
  \item{chr}{integer of length 1 indicating the chromosome number to plot.}
  \item{coord}{integer vector of length 2 containing the start and end coordinates (in
    bases) for the plot.}
  \item{highlight}{(optional) list with two elements: a single numeric value
    \code{coord} and a character \code{strand}. If present, this
    position is marked by a vertical red bar on the coordinate axis.  The color can be changed using the \code{colors} argument below.}
  \item{colors}{(optional) named character vector. If missing, 
    a default color scheme is used:
    \code{c("+"="#00441b", "-"="#081d58", "duplicated"="grey", "cp"="#101010", 
   "highlight"="red", "threshold"="grey")}, 
   where the first three elements refer to the colors of data points and the 
   last three to the colors of lines in the plot.}
  \item{doLegend}{logical: should the plot contain a legend?} 
  \item{featureExclude}{character vector of names of feature types (in
    \code{gff}) that should not be plotted. Default is \code{"chromosome"}, 
    \code{"nucleotide_match"} and \code{"insertion"}.  Additional possibilities 
    include: \code{"ARS"}, \code{"repeat\_region"}, 
    \code{"repeat\_family"} and \code{"nc\_primary\_transcript"}.}
  \item{featureColorScheme}{numeric scalar, used to select a color scheme for the 
   boxes representing genomic features such as coding sequences, ncRNAs etc. 
   Currently the only value supported is 1 (see \code{plotAlongChromLegend} or \code{plotFeatures} for further information).}
  \item{extras}{a matrix containing additional values to be plotted along the 
    chromosome in a separate panel (such as p-values).  This option is only 
    available when \code{y} is specified.  These values should be on the scale [0,1].}
  \item{rowNamesHeatmap}{character vector of row names for the main heatmap.}
  \item{rowNamesExtras}{character vector of row names for the extra heatmap.}
  \item{ylab}{character label for y-axis of main plot.}
  \item{ylabExtras}{character label for y-axis on \code{extras} panel (if specified).}
  \item{main}{character: plot title.}
  \item{colHeatmap}{function describing color scheme for the main heatmap plot (defaults to
    \code{YlGnBu} from RColorBrewer package).}
  \item{colExtras}{function describing color scheme for the extra heatmap plot (if specified)   
    (defaults to \code{Reds} from RColorBrewer package).}
  \item{sepPlots}{logical scalar.  If TRUE, each column of intensities in \code{segObj} or 
    \code{y} is plotted separately (maximum of 3) in the same figure. When FALSE, the average 
    is plotted.  This argument is only used when \code{what} is set to \code{dots}.}
  \item{reOrder}{logical scalar (only used when sepPlots is TRUE). If TRUE, the first 
    column of intensities is printed at the bottom of each plot, and the subsequent columns 
    are plotted above.  If FALSE, the first appears at the top, and the subsequent columns 
    are plotted below.}
  \item{...}{further arguments that can be passed to the
    functions that implement the \code{what} option above (see
    \code{plotSegmentationDots} and 
    \code{plotSegmentationHeatmap}) or \code{gff} plotting (see \code{plotFeatures} and \code{plotAlongChromLegend}).}
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

  \emph{\dots and the different signal visualization methods (\code{what} option)}:
  If \code{what=="dots"}, the argument \code{showConfidenceIntervals}
  can  be a logical scalar to choose whether vertical dashed lines are
  drawn for the confidence interval. In any case, these are only drawn
  if they are present in the \code{segmentation} object in \code{segObj}.
} 
  

\author{Wolfgang Huber <huber@ebi.ac.uk>}
\examples{
  ## 1. see viewSegmentation.R script in the inst/scripts directory
  ## 2. (newer): segmentation.Rnw
  ## 3. (newer): see the plotALongChrom vignette
  data(segnf)
  data(gffSub)
  nmLabel = colnames(segnf$"1.+"@y)
  plotAlongChrom(segnf,chr=1,coord=c(35000,50000),
  			gff=gffSub,rowNamesHeatmap=nmLabel) ##the dots
  plotAlongChrom(segnf,chr=1,coord=c(35000,50000),what="heatmap", 
  			gff=gffSub,rowNamesHeatmap=nmLabel) ##the heatmap
  			
  plotAlongChrom(segnf,chr=1,coord=c(35000,50000),gff=gffSub,
  			showConfidenceIntervals=FALSE) ##do not show the segment confidence interval	
   }

