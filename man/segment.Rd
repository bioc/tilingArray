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
\code{length(x)*maxk*maxcp} in time.}

\value{
  An object of class \code{\link[segmentation-class]{segmentation}}.
}

\author{W. Huber \email{huber@ebi.ac.uk}}

\references{A statistical approach for CGH microarray data analysis.
  Franck Picard, Stephane Robin, Marc Lavielle, Christian Vaisse,
  Gilles Celeux, Jean-Jacques Daudin,
  Rapport de recherche No. 5139, Mars 2004,
  Institut National de Recherche en Informatique et en Automatique (INRIA),
  ISSN 0249-6399.
  The code of this function is based on the Matlab
  implementation presented at  
  \url{http://www.inapg.fr/ens_rech/mathinfo/recherche/mathematique/outil.html},
  but it has evolved.}

\examples{
  x = rep( sin((0:4)/2*pi), each=3) + rnorm(3*5, sd=0.1)
  res = segment(x, maxseg=6, maxk=15)
}
\keyword{manip}
