\name{qualityReport}
\alias{qualityReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Assess quality of normalized poly-A-RNA samples}
\description{
  Function to analyze quality of (normalized) tiling array
  data. Provides plots of intensities along the chromosome around
  specified genes, scatterplots comparing raw and normalized data, and a
  histogram of z-scores for intensities within annotated coding DNA regions.
}
\usage{
qualityReport(x, hybeType, normRef = NULL, compress = TRUE, selectGenes, gff, probeAnno, output = "HTML", outputDir = ".", verbose = TRUE)
}
\arguments{
  \item{x}{tiling array data,
    either a character vector, naming the CEL files to analyze, OR
    an \code{\link[affy]{AffyBatch}} containing the unnormalized data}
  \item{hybeType}{character string, naming the type of hybridized data,
    either 'direct' for directly hybridized RNA (no normalization
    necessary) or 'reverse' for reverse-transcribed RNA}
  \item{normRef}{samples used for normalization of the data, either a
    character vector naming the CEL-files to use for normalization of
    \code{x} of an \code{\link[affy]{AffyBatch}} containing the expression
    data of those samples. \code{NULL} means 'no normalization'}
  \item{compress}{logical; if \code{x} and \code{normRef} were provided as
    character vectors, are the related CEL-files compressed?}
  \item{selectGenes}{character vector; for which genes are
    along-chromosome plots to be prepared}
  \item{gff}{data frame with genome annotation from the GFF file.}
  \item{probeAnno}{environment with probe annotations, see help of
  \code{\link{plotAlongChrom2}} for details}
  \item{output}{character string, type of output to produce, currently
    only 'HTML' is supported.}
  \item{outputDir}{path to storage directory of quality reports}
  \item{verbose}{logical, do you want extended progess output to STDOUT?}
}
%\details{~~ If necessary, more details than the __description__  above ~~}
\value{
  This function is called for producing quality reports in a
  specified format and does not return anything. 
}
\author{Wolfgang Huber}
\seealso{\code{\link[affy]{AffyBatch}},\code{\link{plotAlongChrom2}}}
\examples{}
\keyword{hplot}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
