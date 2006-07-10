\name{segChrom}
\alias{segChrom}

\title{Fit a piecewise constant curve to along chromosome data}
\description{
  This function calls the \code{segment} function for each strand of one or more chromosomes specified by the user.}

}
  
\usage{
segChrom(y, probeAnno, chr=1, nsegs, nrBasesPerSegment = 1500, 
             step=7, confint=FALSE, confintLevel=0.95, save=FALSE, savedir=NULL)
}
\arguments{  
  \item{y}{eSet containing the data to be segmented.}
  \item{probeAnno}{environment with probe annotations. See
    package \code{davidTiling} for an example (?probeAnno).}
  \item{chr}{integer scalar or vector specifying which chromosome(s) to segment.}
  \item{nsegs}{integer vector specifying the maximum number of segments (= 1 + maximum
    number of change points) for each chromosome and strand.  If missing, this value is calculated by dividing the number of bases on a chromosome by \code{nrBasesPerSegment}.}
  \item{nrBasesPerSegment}{integer (length 1): maximum length of a segment.}
  \item{step}{integer, indicating the minimum distance between the consecutive probes. In cases when probes are offset by less than \code{step} bases, the probes are sampled to achieve the desired spacing.}
  \item{confint}{logical scalar.  If \code{TRUE}, confidence intervals for each change-point are calculated.  Otherwise (\code{FALSE}), no confidence intervals are returned.}
  \item{confintLevel}{numeric scalar between 0 and 1 indicating the probability level for the confidence intervals calculated for each change-point.}
  \item{save}{logical scalar. If \code{TRUE} the results of each segmentation are saved to file in the \code{savedir} directory. If \code{FALSE} (default), the results are not saved to file.}
  \item{savedir}{optional character string specifying the filepath where the segmentation results will be saved. Only used when \code{save=TRUE}. Defaults to current working directory.}

}
\details{This function is a wrapper for the \code{segment} function which fits piecewise constant curves to along chromosome data. It calls \code{segment} for each strand of the chromosomes specified by the \code{chr} argument. Refer to the \code{segment} help page for further details.}

\value{
  An environment with containing S4 objects of class \code{"segmentation"} called "1.+", "1.-", etc (depending on \code{chr}), where "+" and "-" indicates the strand and the preceding number refers to the chromosome number. If \code{save=TRUE}, a series of files, "1.+.rda", "1.-.rda", etc (depending on \code{chr}) are saved in the directory specified.}

\author{Matt Ritchie <ritchie@ebi.ac.uk> and Wolfgang Huber <huber@ebi.ac.uk>}

\examples{## library("davidTiling")
## data("davidTiling")
## data("probeAnno")
## isDNA = seq(1:3)
## yn = normalizeByReference(davidTiling[,-isDNA],davidTiling[,isDNA], probeAnno=probeAnno)
## seg1 = segChrom(yn, probeAnno, chr=1) # this will take a while to run!
}
\keyword{manip}
