\name{findSegments}
\alias{findSegments}
\alias{segment}
\title{Fit a piecewise constant curve to a sequence of numbers  --
  OBSOLETE, please use function segment instead.}
\description{
  This function is only here for backward compatibility - please use
  \code{link{segment}}.
  
  The function fits a piecewise constant curve to a sequence of numbers
  using a simple least squares cost function and the dynamic programming
  algorithm of Picard et al. (full citation see below).
}
\usage{
findSegments(x, maxcp, maxk, verbose=TRUE)
}
\arguments{
  \item{x}{Numeric (real) vector.}
  \item{maxcp}{Integer (length 1): maximum number of segments (= 1 + maximum
    number of change points).}
  \item{maxk}{Integer (length 1): maximum length of a segment.}
  \item{verbose}{Logical: if this parameter has a positive
    value, various diagnostic output is printed.}
}
\details{The complexity of the algorithm is 
\code{length(x)*maxk} in memory and
\code{length(x)*maxk*maxcp} in time.}

\value{
  An object of class \code{"segmentation"}
  A list with elements
  \item{J}{likelihood criterion}
  \item{th}{matrix of segment start points}
  \item{dat}{the data used for the segmentation}
  \item{call}{the function call}.

  See the vignette, and
  the paper cited below for details.
}

\note{This function is depracated and replaced by function
  \code{segment}, but still included for backward compability. 
}

\author{W. Huber \email{huber@ebi.ac.uk}, Joern Toedling
  \email{toedling@ebi.ac.uk}}

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
 res = findSegments(x, maxcp=6, maxk=15)
}
\keyword{manip}
