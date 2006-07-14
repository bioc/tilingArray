PMindex <- function(probeAnno) {
#  as.integer(probeAnno$matchInfo[,"PMindex"])
  nprobes <- length(probeAnno$probeReverse$no_feature)
  isPM = logical(nprobes)
  for (j in probeAnno$probeReverse) isPM[as.character(j) != ""] = TRUE
  seq(1:nprobes)[isPM]
  }

MMindex <- function(probeAnno) {
  as.integer(probeAnno$matchInfo[,"MMindex"])
  }

BGindex <- function(probeAnno) {
  nprobes <- length(probeAnno$probeReverse$no_feature)
  seq(1:nprobes)[(probeAnno$probeReverse$no_feature == "no"
   & probeAnno$probeDirect$no_feature == "no")]
 }
