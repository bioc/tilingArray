\name{scoreSegments}
\alias{scoreSegments}
\title{BLA}
\description{BLA}
\usage{
scoreSegments(x, gff, 
  nrBasePerSeg = 1500, 
  probeLength  = 25,
  knownFeatures = c("CDS", "gene", "ncRNA", "nc_primary_transcript",
        "rRNA", "snRNA", "snoRNA", "tRNA",
        "transposable_element", "transposable_element_gene"),
  verbose = TRUE)
}
\arguments{
  \item{x}{character vector.}
  \item{field}{character vector of length 1, containing the field name.}
  \item{attrsep}{character vector of length 1, containing the separator name.}
}
\details{
## This algorithm scores segments.
## For each segment, we calculate and/or record its:
## * chromosome, strand, start position, end position (in bases)
## * level (mean-background)
## * p-value (one-sample t-test againt hypothesis mean=0)
##
## If it overlaps with an annotated feature
## * the name of that feature
## * the signed distance between its start and that of the feature
## * the signed distance between its end and that of the feature
## If not:
## * the distance between its start and the end of the next feature to the left
## * the distance between its end and the start of the next feature to the left
## 
## and the same information again for the opposite strand
##--------------------------------------------------------------------------------

## For the definition of pseudogenes at SGD, see Docs/PseudogenesAtSGD.pdf
}

\value{
  Character vector.
}

\author{W. Huber <huber@ebi.ac.uk>}

\examples{
}
\keyword{manip}
