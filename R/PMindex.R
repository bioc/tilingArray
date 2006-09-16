PMindex <- function(probeAnno) {
  isPM = logical(length(probeAnno$probeReverse$no_feature))
  for (j in probeAnno$probeReverse)
    isPM[as.character(j) != ""] = TRUE
  which(isPM)
}

BGindex <- function(probeAnno) {
  which(probeAnno$probeReverse$no_feature == "no" & probeAnno$probeDirect$no_feature == "no")
}
