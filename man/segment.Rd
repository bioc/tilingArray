\name{segment}
\alias{segment}

\title{Fit a piecewise constant curve: segmentation by dynamic
  programming}

\description{
  The function fits a piecewise constant curve to one or multiple
  sequences of measurements,
  using a least squares cost function and an O(n) dynamic programming
  algorithm (see references).
}
\usage{
segment(y, maxseg, maxk)
}
\arguments{
  \item{y}{Numeric matrix. Rows correspond to the
      \code{x}-variable, columns to replicate measurements at the same
      value of \code{x}. Breakpoints are fitted along the \code{x}-axis.
      For example, the \code{x}-variable can be genomic coordinates or
      time. The segmentation will be along the rows of \code{y}.}
  \item{maxseg}{integer of length 1, maximum number of segments (= 1 + maximum
    number of change points).}
  \item{maxk}{integer of length 1, maximum length of a single segment.}
}
\details{The complexity of the algorithm is 
\code{length(x)*maxk} in memory and
\code{length(x)*maxk*maxseg} in time.}

\value{
  An object of class \code{\link[tilingArray:class-segmentation]{segmentation}}.
}

\author{W. Huber \email{huber@ebi.ac.uk}}

\references{
  [1] Transcript mapping with high-density oligonucleotide tiling
  arrays. Huber W, Toedling J, Steinmetz, L.
  Bioinformatics 22, 1963-1970 (2006).
  
  [2] A statistical approach for CGH microarray data analysis.
  Franck Picard, Stephane Robin, Marc Lavielle, Christian Vaisse,
  Gilles Celeux, Jean-Jacques Daudin.
  BMC Bioinformatics. 2005 Feb 11; 6:27.
}

\examples{
  x = rep( sin((0:4)/2*pi), each=3) + rnorm(3*5, sd=0.1)
  res = segment(x, maxseg=6, maxk=15)
}
\keyword{manip}
