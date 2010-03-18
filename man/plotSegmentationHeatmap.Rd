\name{plotSegmentationHeatmap}
\alias{plotSegmentationHeatmap}
\title{Plot a heatmap diagram for a region along a chromosome}
\keyword{hplot}
\description{Plot a heatmap diagram for a region along a chromosome}

\usage{
plotSegmentationHeatmap(dat, xlim, ylab, rowNames, 
                        chr=1, strand="+", vpr, colors, 		
                        colHeatmap=colorRamp(brewer.pal(9, "YlGnBu")), 
                        showConfidenceIntervals=TRUE,
                        just=c("left","centre"),
                        main,makeRasterImage = TRUE, ...)
}

\arguments{
  \item{dat}{list containing data to be plotted (see \emph{Details} section below for particulars).}
  \item{xlim}{integer vector of length 2 with start and end coordinates (in bases) for plotting.}
  \item{ylab}{character scalar specifying y-axis label.}
  \item{rowNames}{character vector specifying a name for each row in the heatmap plot.}
  \item{chr}{integer of length 1 indicating the chromosome to plot (defaults to 1).}
  \item{strand}{character scalar which should be set to either \code{+} or \code{-} to 
    indicate which strand of DNA to plot the intensity values from (defaults to "+").}
  \item{vpr}{which viewport to plot the figure in.  If this function is called directly 
   by the user this argument should be left missing.}
  \item{colors}{named character vector, optional. If missing, 
    a default color scheme is used:
    \code{c("+"="#00441b", "-"="#081d58", "duplicated"="grey", "cp"="#101010", 
   "highlight"="red", "threshold"="grey")}, 
   where the first three elements refer to colors of data points and the 
   last three to those of lines in the plot.}
  \item{colHeatmap}{function describing color scheme for the heatmap plot (defaults to
    \code{YlGnBu} from RColorBrewer package).}  
  \item{showConfidenceIntervals}{logical scalar indicating whether
    confidence intervals for each change-point are to be plotted (only
    available once segmentation has occurred).}
  \item{just}{character vector specifying the justification of the
    supplied values to the given coordinates; setting the first entry to
    "left" indicates that the supplied x-coordinates are the start
    positions of the probes, change this to "centre" if the
    x-coordinates are the probe middle positions. Usually the second
    entry should be "centre" (see \code{\link[grid]{grid.rect}})}
  \item{main}{character vector specifying plot title.}
  \item{makeRasterImage}{logical scalar indicating whether to plot the heatmap image
  	by the grid.raster (see \code{\link[grid]{grid.raster}}) or 
  	the grid.rect (see \code{\link[grid]{grid.rect}}) function in grid package.The default is to 
  	generate raster image, as it can be displayed much faster with a relatively smaller file size.}
  \item{...}{additional arguments.}
}

\details{
  This function is called by \code{plotAlongChrom} if the argument
  \code{what} is set to \code{heatmap}.
  Although this function can be called directly by the user, this is not recommended.
  The \code{dat} list contains the following items:
  \describe{
    \item{\code{x}}{x-coordinates (in bases) along chromosome}
    \item{\code{y}}{intensity matrix of probes along chromosome}
    \item{\code{flag}}{indicates probe uniqueness in the genome.  Possibilities are 3: multiple perfect matches, 2: has no PM but one or more near-matches, 1: has exactly one PM and some near-matches in the genome, 0: has exactly one PM and no near-matches.}
    \item{\code{extras}}{(optional) matrix of additional values (such as test-statistics/p-values) to be plotted}
  }
}
\examples{
  data(segnf)
  data(gffSub)
  nmLabel = colnames(segnf$"1.+"@y)
  plotAlongChrom(segnf,chr=1,coord=c(35000,50000),what="heatmap", 
  		gff=gffSub,rowNamesHeatmap=nmLabel) ##using raster image
}
\author{Wolfgang Huber <huber@ebi.ac.uk>}
