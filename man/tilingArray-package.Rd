\name{tilingArray-package}
\alias{tilingArray-package}
\alias{tilingArray}
\docType{package}
\title{tilingArray package overview}
\description{
   tilingArray package overview
}
\details{
  The package provides some functionalities that can be useful for the
  analysis of high-density tiling microarray data (such as Affymetrix
  genechips) for measuring transcript abundance and architecture.
  The main functionalities of the package are:
  
  \itemize{
    \item The segmentation class
     for representing partitionings of a linear series of data (such as
     microarray intensity readings along a chromosome strand).
   \item The function \code{\link{segment}}
     for fitting piecewise constant models using a dynamic programming
     algorithm that is both fast and exact,
     and \code{\link[segmentation-class]{confint}} for
     calculating confidence intervals using the
     strucchange package.
     Please see the vignette \emph{Segmentation demo} in the file
     inst/doc/segmentation.pdf (source file inst/scripts/segmentation.Rnw).
   \item The function \code{\link{plotAlongChrom}} for generating pretty
     plots of segmentations along with genomic features. Please also see
     the vignette \emph{Segmentation demo}.
   \item The function \code{\link{normalizeByReference}}
     for probe-sequence dependent response adjustment from 
     a (set of) reference hybridizations.
     Please see the vignette \emph{Assessing signal/noise ratio before and after normalization} in the file
     inst/doc/assessNorm.pdf (source file inst/scripts/assessNorm.Rnw).
     }
}
\author{
   W. Huber <huber@ebi.ac.uk>
}
\keyword{package}
