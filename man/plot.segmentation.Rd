\name{plot.segmentation}
\alias{plot.segmentation}
\title{Plot likelihood versus model complexity penalty}
\keyword{hplot}
\description{Plot likelihood versus model complexity penalty.}

\usage{plot.segmentation(x, ...)}

\arguments{
  \item{x}{object of S3 class \code{segmentation}}
  \item{...}{arguments that are passed to \code{plot.default}}
}

 \details{No details.}

\author{Wolfgang Huber <huber@ebi.ac.uk>}
\examples{
 x = rep( sin((0:4)/2*pi), each=3) + rnorm(3*5, sd=0.1)
 res = findsegments(x, maxcp=6, maxk=20)
 if(interactive())
  plot(res)
 }

