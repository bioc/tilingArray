\name{plotAlongChromLegend}
\alias{plotAlongChromLegend}
\title{Plot a legend for genomic features}
\keyword{hplot}
\description{Plot a legend for genomic features}

\usage{
plotAlongChromLegend(vpr, nr=2, 
    featureColorScheme=1, 
    featureExclude=c("chromosome", "nucleotide_match", "insertion"),
    mainLegend, cexLegend=0.35, cexMain=1)
}

\arguments{
  \item{vpr}{vector specifying where to place the legend in figure (set up by using the 
   \code{viewport} function from the grid package.  When this function is called directly 
   by the user this argument should be left missing.}
  \item{nr}{numeric scalar, specifying the number of rows to plot legend over 
   (default value is 2).}
  \item{featureColorScheme}{numeric scalar, used to select a color scheme for the 
   boxes representing genomic features such as coding sequences, ncRNAs etc. 
   Currently the only value supported is 1.}
  \item{featureExclude}{character vector of names of feature types (in
   \code{gff}) that should not be plotted. Default is \code{"chromosome"}, \code{"nucleotide_match"} and 
   \code{"insertion"}.  Additional possible candidates include: \code{"ARS"}, 
   \code{"repeat\_region"}, \code{"repeat\_family"} and 
   \code{"nc\_primary\_transcript"}.}
  \item{mainLegend}{character vector specifying legend title.}
  \item{cexLegend}{numeric scalar specifying the magnification to be used for the legend 
   text relative to the current text size.}
  \item{cexMain}{numeric scalar specifying the magnification to be used for the legend 
   title relative to the current text size.}
}

\details{
  This function is usually called by \code{plotAlongChrom} when \code{doLegend} is TRUE.  It can also be called directly by the user to produce a separate legend.  

The following features are included in the legend (unless excluded using the \code{featuredExclude} option): \code{"chromosome"}, \code{"nucleotide_match"}, \code{"pseudogene"}, \code{"uORF"}, \code{"nc_primary_transcript"}, \code{"region"}, \code{"repeat_family"}, \code{"repeat_region"}, \code{"transposable_element"}, \code{"transposable_element_gene"}, \code{"ARS"}, \code{"centromere"}, \code{"telomere"}, \code{"insertion"}, \code{"CDS"}, \code{"CDS_dubious"}, \code{"ncRNA"}, \code{"tRNA"}, \code{"snRNA"}, \code{"rRNA"}, \code{"snoRNA"}, \code{"binding_site"} and \code{"TF_binding_site"}.}

\author{Wolfgang Huber <huber@ebi.ac.uk>}
\examples{## plotAlongChromLegend(mainLegend="Legend")}
