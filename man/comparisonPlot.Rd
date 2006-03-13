\name{comparisonPlot}
\alias{comparisonPlot}

\title{
  Plot a vertical layout of panels for the comparison of different
  along-chromosome profiles.
}
\description{
  This function is used for Figure 5 in the David et al. (PNAS 2006)
  paper and in the Huber et al. () paper.
}
\usage{
comparisonPlot(x, y, xscale=range(x), yscale, anno, ticks)
}
\arguments{
  \item{x}{numeric vector.}
  \item{y}{list of numeric vector, each of same length as \code{x}.}
  \item{xscale}{numeric vector of length 2.}
  \item{yscale}{matrix with 2 rows and columns corresponding to the
    elements of \code{x}.}
  \item{anno}{dataframe with columns \code{start}, \code{end}, and
    \code{name}, each row corresponds to one gene CDS to be plotted at the bottom.}
  \item{ticks}{numeric vector, where to plot the ticks.}
}
\details{}

\value{
  Function is called for its side-effect.
}

\author{W. Huber \email{huber@ebi.ac.uk}}

\references{...}

\examples{
  ##
}
\keyword{manip}