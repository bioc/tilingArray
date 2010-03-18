\name{gffSub}
\alias{gffSub}
\title{Example of a genomic feature object }
\description{Example of a genomic feature object}
\usage{
 data(gffSub)
}
\details{
  gffSub is a data frame that contains the features from 35000bp - 50000bp in yeast chromosome one.    	  
}

\author{Zhenyu Xu <zhenyu@ebi.ac.uk>}

\examples{
  data(segnf)
  data(gffSub)
  nmLabel = colnames(segnf$"1.+"@y)
  plotAlongChrom(segnf,chr=1, coord=c(35000,50000),what="heatmap", 
  		gff=gffSub,rowNamesHeatmap=nmLabel)
}
\keyword{data}


