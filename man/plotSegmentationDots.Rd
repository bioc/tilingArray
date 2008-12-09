\name{plotSegmentationDots}
\alias{plotSegmentationDots}
\title{Plot points for a region along a chromosome}
\keyword{hplot}
\description{Plot points for a region along a chromosome}

\usage{
plotSegmentationDots(dat, xlim, ylim, ylab, threshold=NA, 
                     chr=1, strand="+", vpr, colors, main, 
                     pointSize=unit(0.6, "mm"), showConfidenceIntervals=TRUE, 
                     sepPlots=FALSE, cexAxisLabel=1, cexAxis=1,...)
}

\arguments{
  \item{dat}{list containing data to be plotted (see \emph{Details} section below for particulars).}
  \item{xlim}{integer vector of length 2 with start and end coordinates (in bases) for plotting.}
  \item{ylim}{numeric vector containing the y limits of the plot.}
  \item{ylab}{character scalar (if \code{sepPlots=FALSE}) or vector containing y-axis label(s).}
  \item{threshold}{numeric scalar indicating the threshold of expression (default value is NA, for 
   no threshold.  If a value is supplied, it is subtracted from the intensity measures in \code{dat$y}.}
  \item{chr}{integer of length 1 indicating the chromosome to be plot (defaults to 1).}
  \item{strand}{character scalar which should be set to either \code{+} or \code{-} to 
    indicate which strand of DNA to plot the intensity values from (defaults to "+").}
  \item{vpr}{which viewport to plot the figure in.  If this function is called directly 
   by the user this argument should be left missing.}
  \item{colors}{named character vector, optional. If missing, 
    a default color scheme is used:
    \code{c("+"="#00441b", "-"="#081d58", "duplicated"="grey", "cp"="#101010", 
   "highlight"="red", "threshold"="grey")}, 
   where the first three elements refer to the colors of data points and the 
   last three to the colors of lines in the plot.}
  \item{main}{character vector specifying plot title.}
  \item{pointSize}{an object of class unit which specifies the size of each point.  Default value is \code{unit(0.6, "mm")}.}
  \item{showConfidenceIntervals}{logical scalar indicating whether confidence intervals for each change-point are to be plotted (only available once segmentation has occurred).}
  \item{sepPlots}{logical scalar indicating whether the intensities are plotted separately for each array (if \code{dat$y} has multiple columns).  Defaults to FALSE, in which case the average intensity for each probe is plotted.  When TRUE, up to 3 arrays can be plotted separately (more than 3 gets crowded).}
  \item{cexAxisLabel}{numeric scalar specifying the magnification to be used for the y-axis label relative to the current test size.}
  \item{cexAxis}{numeric scalar specifying the magnification to be used for the y-axis annotation relative to the current text size.}
  \item{...}{additional arguments.}
}

\details{
  This function is called by \code{plotAlongChrom} when the argument \code{what} is set to \code{dots}.  Although this function can be called directly by the user, this is not recommended.

The \code{dat} list contains the following items:
items \code{x}: x-coordinates (in bases) along chromosome, 
\enumerate{
    \item \code{y}: intensity matrix of probes along chromosome,
    \item \code{flag}: indicates probe uniqueness in the genome.  Possibilities are 3: multiple perfect matches, 2: has no PM but one or more near-matches, 1: has exactly one PM and some near-matches in the genome, 0: has exactly one PM and no near-matches.
    \item \code{extras}: (optional) matrix of additional values (such as test-statistics/p-values) to be plotted.
    }
}
      
\author{Wolfgang Huber <huber@ebi.ac.uk>}
