\name{findsegments}
\alias{findsegments}
\title{Fit a piecewise constant curve to a sequence of numbers}
\description{
  The function fits a piecewise constant curve to a sequence of numbers
  using a simple least squares cost function and the dynamic programming
  algorithm of Picard et al. (full citation see below).
}
\usage{
findsegments(x, maxcp, maxk, verbose=0)
}
\arguments{
  \item{x}{Numeric (real) vector.}
  \item{maxcp}{Integer (length 1): maximum number of segments (= 1 + maximum
    number of change points).}
  \item{maxk}{Integer (length 1): maximum length of a segment.}
  \item{verbose}{Integer (length 1): if this parameter has a positive
    value, various diagnostic output is printed.}
}
\details{The complexity of the algorithm is 
  \code{length(x) * maxk}.
}

\value{
  A list with two elements \code{J} and \code{th}. See the below
  cited paper for details (and soon it will also be explained
  in the vignette, please bear with me or remind me).
}

\author{W. Huber <huber@ebi.ac.uk>}

\references{A statistical approach for CGH microarray data analysis.
  Franck Picard, Stephane Robin, Marc Lavielle, Christian Vaisse,
  Gilles Celeux, Jean-Jacques Daudin,
  Rapport de recherche No. 5139, Mars 2004,
  Institut National de Recherche en Informatique et en Automatique (INRIA),
  ISSN 0249-6399.
  The code of this function is based on the Matlab
  implementation presented at  
  http://www.inapg.fr/ens_rech/mathinfo/recherche/mathematique/outil.html,
  but it has evolved.}

\examples{
 x = rep( sin((0:4)/2*pi), each=3) + rnorm(3*5, sd=0.1)
 res = findsegments(x, maxcp=6, maxk=20)
}
\keyword{manip}
