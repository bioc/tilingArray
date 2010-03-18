\name{segnf}
\alias{segnf}
\title{Example of a segmentation output object }
\description{Example of a segmentation output object}
\usage{
  data(segnf)
}
\details{
  segnf is an environment which contains two segmentation object for the "+"
  and "-" strand of yeast chromosome one. Each segmentation object contains
  the tiling array expression profiling between condition YPE and YPD, 3 replicates each.
  The raw array data are normalized, mapped to the yeast genome and segmented using the
  segChrom function. The object is the output of segChom with some small modification.
  Only a small region (35000bp-50000bp) is included in this data     	  
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


