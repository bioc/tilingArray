\name{plotFeatures}
\alias{plotFeatures}
\title{Plot genomic features for a region along a chromosome}
\keyword{hplot}
\description{Plot genomic features for a region along a chromosome}

\usage{
plotFeatures(gff, chr, xlim, strand, vpr, featureColorScheme=1, 
             featureExclude=c("chromosome", "nucleotide_match", "insertion"),            
             featureNoLabel=c("uORF", "CDS"), ...)
}

\arguments{
  \item{gff}{data frame with genome annotation from the GFF file.}
  \item{chr}{integer of length 1 specifying the chromosome to plot the features for.}
  \item{xlim}{integer of length 2 with start and end coordinates (in bases) for plotting.}
  \item{strand}{character scalar which should be set to either \code{+} or \code{-} to 
    indicate which strand of DNA to plot the features from.}
  \item{vpr}{which viewport to plot the features in.}
  \item{featureColorScheme}{numeric scalar, used to select a color scheme for the 
   boxes representing genomic features such as coding sequences, ncRNAs etc. 
   Currently the only value supported is 1.}
  \item{featureExclude}{character vector of names of feature types (in
    gff) that should not be plotted. Default is \code{"chromosome"}, \code{"nucleotide_match"} and 
    \code{"insertion"}.  Additional possible candidates include: \code{"ARS"}, \code{"repeat\_region"},
    \code{"repeat\_family"} and \code{"nc\_primary\_transcript"}.}
  \item{featureNoLabel}{character vector, names of feature types (in
    gff) that should not be labelled with their names (if they are plotted).}
  \item{...}{additional arguments.}
}

\details{
  This function is called by \code{plotAlongChrom} when the \code{gff} argument has been specified.  It should not be called directly by the user.}
      
\author{Wolfgang Huber <huber@ebi.ac.uk>}
\examples{}
