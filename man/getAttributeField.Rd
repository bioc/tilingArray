\name{getAttributeField}
\alias{getAttributeField}
\title{Extract the value of a certain field out of a character vector
  such as in the "attributes" column of a GFF table.}
\description{
  Extract the value of a certain field out of a character vector
  such as in the "attributes" column of a GFF table.
}
\usage{
getAttributeField(x, field, attrsep=";") 
}
\arguments{
  \item{x}{character vector.}
  \item{field}{character vector of length 1, containing the field name.}
  \item{attrsep}{character vector of length 1, containing the separator name.}
}
\details{See example.
}

\value{
  Character vector.
}

\author{W. Huber <huber@ebi.ac.uk>}

\examples{
  acol = c("ID=46891;Name=TEL01L-TR;Note=Bla",
   "ID=46892;Name=TEL01L;Note=Di",
   "ID=46893;Name=TEL01L-XR;Note=Bla")

  getAttributeField(acol, "Name")
  getAttributeField(acol, "ID")
}
\keyword{manip}
