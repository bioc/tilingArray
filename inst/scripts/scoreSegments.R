options(error=recover, warn=2)
library("tilingArray")
source("/homes/huber/madman/Rpacks/tilingArray/R/scoreSegments.R")

rnaTypes  = c("seg-polyA-050525", "seg-tot-050525", "seg-tot2-050525")[3]
doNotLoadSegScore=TRUE
source("scripts/readSegments.R")

##if(!exists("xn")) {
##  fn = "seg-dir-050521/xn.rda"
##  cat("Loading", fn, "\n")
##  load(fn)
##}

addDirectHybe = function(s) {
  dhl = numeric(nrow(s))
  strand = otherStrand(s$strand)
  y = exprs(xn)[, "050507_dirRNA_10ug_F1.cel.gz"]
  cat("Adding direct hybe: ")
  for(i in 1:nrow(s)) {
    if(i%%1000==0)cat(i, "")
    uni = get(paste(s$chr[i], strand[i], "unique", sep="."), probeAnno)
    ind = get(paste(s$chr[i], strand[i], "index",  sep="."), probeAnno)
    sta = get(paste(s$chr[i], strand[i], "start",  sep="."), probeAnno)
    end = get(paste(s$chr[i], strand[i], "end",    sep="."), probeAnno)
    sel = ind[uni & (sta >= s[i, "start"]) & (end <= s[i, "end"])]
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
    save(segScore, file=file.path(indir[rt], sprintf("segScore-%d.rda", as.integer(nrbps))),
         compress=TRUE)
  } 
}
