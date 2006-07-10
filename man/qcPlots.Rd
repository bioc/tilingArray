\name{qcPlots}
\alias{qcPlots}
\title{Generate simple diagnostic plots for Affymetrix tiling array data}
\keyword{hplot}
\description{Generate simple diagnostic plots for Affymetrix tiling array data}

\usage{
qcPlots(x, html=TRUE, plotdir=NULL, probeAnno, gff, 
                    chr=4, coord=c(230000,245000),
  		    nr = 2560, nc = 2560,
		    ylimchrom=c(5,16), nucleicAcid, pmindex, pgm=TRUE,
                    ext=".cel", ranks=FALSE, ...)
}

\arguments{
  \item{x}{eSet containing the data to be plotted.}
  \item{html}{logical scalar.  If \code{TRUE} an html summary page 
              'qcsummary.htm' is generated.  If \code{FALSE}, no summary page is generated.}
  \item{plotdir}{optional character string specifying the filepath where the plots will be saved.  Defaults to current working directory.}
  \item{probeAnno}{environment with probe annotations. See
    package \code{davidTiling} for an example (?probeAnno).}
  \item{gff}{data frame with genome annotation from the GFF file.} 
  \item{chr}{integer of length 1 indicating the chromosome number to plot.}
  \item{coord}{integer vector of length 2 containing the start and end coordinates (in
    bases) for the along chromosome intensity plot.}
  \item{nr}{integer, indicating the number of probes in each row on the array (2560 for yeast tiling arrays).}
  \item{nc}{integer, indicating the number of probes in each column on the array (2560 for yeast tiling arrays).}
  \item{ylimchrom}{numeric vector containing the y limits of the along chromosome intensity plot.}
  \item{nucleicAcid}{character vector or factor indicating what sample has been hybridised to each
array. Used to color the boxplots and smoothed histograms of intensities.}
  \item{pmindex}{integer vector of indices of PM probes in \code{x}.  If missing, this information is extracted from \code{probeAnno}.}
  \item{pgm}{logical scalar.  If \code{TRUE}, image plots will be saved as .pgm files.  Otherwise (\code{FALSE}), they are converted to jpegs. On windows machines, this argument should be set to \code{TRUE}.}
  \item{ext}{character string indicating the file extension.}
  \item{ranks}{logical scalar.  If \code{TRUE}, imageplots will show ranks of standardised probe intensities.  Otherwise (\code{FALSE}, default), the standardised probe intensities are plotted.}
  \item{...}{further arguments that can be passed to the plotting function \code{plotSegmentationDots}.}
}

\details{
  This function creates boxplots, smoothed histogram (density) plots, imageplots and 
  along chromosome plots of the raw (log base 2) probe intensity data.
  
  An html page called 'qcsummary.htm' which displays the results, is created when \code{html=TRUE}. 

  Imageplots of standardised intensities (i.e. (probe intensity - minimum probe intensity) divided by the difference between the maximum and minimum probe intensities, all on log base 2 scale) or the ranks of these standardised intensities are plotted depending on the \code{ranks} argument. 

  The individual plots are named by replacing the file extension (specified by \code{ext}) of each 'celfile.ext', with 'density.png' for smoothed histogram plots, 'gencoord.jpg', for along chromosome plots and either 'log.pgm' ('log.jpg' if \code{pgm=FALSE}) or 'rank.pgm' ('rank.jpg' if \code{pgm=FALSE}) for the imageplots, depending on the \code{ranks} argument.
}

\author{Matt Ritchie <ritchie@ebi.ac.uk> and Wolfgang Huber <huber@ebi.ac.uk>}
\examples{## library(davidTiling)
## data(davidTiling)
## data(probeAnno)
## qcPlots(davidTiling, probeAnno)
 }

