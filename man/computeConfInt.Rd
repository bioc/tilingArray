\name{computeConfInt}
\alias{computeConfInt}
\title{Compute Confidence Intervals for Breakpoints}
\description{
  Computes confidence intervals for breakpoints.
  Does not need to be called by the \code{tilingArray} user.
}
\usage{
    computeConfInt( object, parm = NULL, level = 0.95,
    breaks = NULL, het.reg = FALSE, het.err = TRUE, ...)
}
\arguments{
  \item{object}{a list containing the data, formula, breakpoints, and
    residuals}
  \item{parm}{the same as \code{breaks}, only one of the two should be
    specified.}
  \item{level}{the confidence level required.}
  \item{breaks}{an integer specifying the number of breaks to be used.
    By default the breaks of the minimum BIC partition are used.}
  \item{het.reg}{logical. Should heterogenous regressors be assumed? If set
    to \code{FALSE} the distribution of the regressors is assumed to be
    homogenous over the segments.}
  \item{het.err}{logical. Should heterogenous errors be assumed? If set
    to \code{FALSE} the distribution of the errors is assumed to be
    homogenous over the segments.}
  \item{\dots}{\emph{currently not used}.}
}

\details{
  As the breakpoints are integers (observation numbers) the corresponding
  confidence intervals are also rounded to integers.

  The distribution function used for the computation of confidence
  intervals of breakpoints is given in Bai (1997). The procedure, in
  particular the usage of heterogenous regressors and/or errors, is
  described in more detail in Bai & Perron (2003).

  The breakpoints should be computed from a formula with \code{breakpoints},
  then the confidence intervals for the breakpoints can be derived by
  \code{confint} and these can be visualized with the segmentation.
  For an example see \code{\link{plot.segmentation}}.
}

\value{
  A matrix containing the breakpoints and their lower and upper
  confidence boundary for the given level.
}

\author{Achim Zeileis \email{Achim.Zeileis@R-project.org},
  with minor contributions from Joern Toedling
  \email{toedling@ebi.ac.uk}
  and Wolfgang Huber \email{huber@ebi.ac.uk}}

\references{
Bai J. (1997), Estimation of a Change Point in Multiple Regression Models,
  \emph{Review of Economics and Statistics}, \bold{79}, 551-563.

  Bai J., Perron P. (2003), Computation and Analysis of Multiple
  Structural Change Models,
  \emph{Journal of Applied Econometrics}, \bold{18}, 1-22.
}

\seealso{\code{\link{findSegments}},
         \code{\link[strucchange]{confint.breakpointsfull}}}
       
\examples{
    ### see 'plot.segmentation' for an example with plot
}

\keyword{internal}
