options(error=recover, warn=2)
library("tilingArray")
source("/homes/huber/madman/Rpacks/tilingArray/R/scoreSegments.R")

doNotLoadSegScore=TRUE
source("scripts/readSegments.R")

if(!exists("x")) {
  cat("Loading x.rda\n")
  load("x.rda")
}

addDirectHybe = function(s) {
  dhl = numeric(nrow(s))
  strand = otherStrand(s$strand)
  y = exprs(x)[, "050507_dirRNA_10ug_F1.cel.gz"]
  cat("Adding direct hybe: ")
  for(i in 1:nrow(s)) {
    if(i%%1000==0)cat(i, "")
    ind = get(paste(s$chr[i], strand[i], "index", sep="."), probeAnno)
    sta = get(paste(s$chr[i], strand[i], "start", sep="."), probeAnno)
    sel = ind[(sta >= s[i, "start"]) & (sta <= s[i, "end"])]
    dhl[i] = mean(y[sel])
  }
  cat("\n")
  s$directHybeLevel = dhl
  return(s)
}

nrbpsList = c(1500)

for(rt in rnaTypes) {
  for(nrbps in nrbpsList) {
    cat(">> ", rt, nrbps, "<<\n")
    segScore = scoreSegments(get(rt), gff=gff, nrBasePerSeg=nrbps)
    ## segScore = addDirectHybe(segScore)
    cat("NOT ADDING DIRECT HYBE\n")
    save(segScore, file=file.path(indir[rt], sprintf("segScore-%d.rda", as.integer(nrbps))),
         compress=TRUE)
  } 
}
