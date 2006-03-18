\name{readCel2eSet}
\alias{readCel2eSet}
\title{Read celfiles into an eSet object.}
\description{
  This is a wrapper for \code{\link[affy:ReadAffy]{ReadAffy}} that
  returns an \code{\link[Biobase:eSet]{eSet}} object rather than an
  AffyBatch. This is particularly usefiles for arrays for which we have
  or need no CDF environment.
}
\usage{
readCel2eSet(filenames, path=".", rotated=FALSE)
}
\arguments{
  \item{filenames}{Character vector with CEL file names.}
  \item{path}{Character scalar with path to CEL files.}
  \item{rotated}{Logical scalar, see details.}
}
\details{
  The \code{rotate} options allows to deal with different versions of
  the scanner software. Older versions rotated the image by 90 degrees,
  newer ones do not. Use the default \code{rotated=FALSE} for CEL files
  produced by the newer version.
}

\value{
  Matrix with \code{maxk} rows and \code{length(x)} columns.
}

\author{W. Huber}

\examples{
## To test the rotation, look at the scatterplot between two DNA hybes
## that were measured with scanner software that rotated (041120) and did
## not rotate (060125)
##
## cp /ebi/research/huber/Projects/tilingArray/Celfiles/041120_S96genDNA_re-hybe.cel.gz ~/p/tmp
## cp /ebi/research/huber/Projects/allelicTranscription/celfiles_allelictrans/060125_S96_genomicDNA.zip ~/p/tmp
## cd ~/p/tmp
## gunzip 041120_S96genDNA_re-hybe.cel.gz
## unzip 060125_S96_genomicDNA.zip
##
\dontrun{
library("affy")
library("geneplotter")
options(error=recover)

e1 = readCel2eSet("041120_S96genDNA_re-hybe.cel", rotated=TRUE)
e2 = readCel2eSet("060125_S96_genomicDNA.CEL")

smoothScatter(log(exprs(e1)), log(exprs(e2)), nrpoints=0)
}
}
\keyword{manip}