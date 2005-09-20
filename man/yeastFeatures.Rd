\name{yeastFeatures}
\alias{yeastFeatures}
\non_function{}
\title{A data.frame with metadata about the features in SGD GFF files}
\description{The rows of the data frame correspond to feature
  categories, such as "gene", "CDS", "telomere". The 
  column \code{isTranscribed} is a logical vector that denotes whether
  this feature category is considered to be transcribable.}
\usage{
data(yeastFeatures)
}
\format{A data.frame}
\author{W. Huber <huber@ebi.ac.uk>}
\examples{
data(yeastFeatures)
}
\keyword{datasets}
