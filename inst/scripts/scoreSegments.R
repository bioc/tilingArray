options(error=recover, warn=2)
library("tilingArray")
source("/homes/huber/madman/Rpacks/tilingArray/R/scoreSegments.R")

if(!exists("gff"))
  load("probeAnno.rda")

chrs = 1:17
indirList = c("segmentation-3polyA", "seg-polyA-050428", "seg-tot-050421")[2:3]
nrbpsList = c(1500, 2000)[1]

for(indir in indirList) {
  cat(indir, "\n")

  if(TRUE) {
    s  = new.env()
    cat("Loading ")
    for(chr in chrs) {
      cat(chr, "")
      for(strand in c("+", "-")) {
        fn = file.path(indir, paste(chr, strand, "rda", sep="."))
        load(fn)
        assign(paste(chr, strand, "seg", sep="."), seg, envir=s)
        assign(paste(chr, strand, "dat", sep="."), dat, envir=s)
      }
    } ## for chr
    cat("\n")
  } else {
    cat("NOT LOADING DATA FILES!\n")
  }
  
  for(nrbps in nrbpsList) {
    cat(">> ", nrbps, "<<\n")
    segScore = scoreSegments(s, gff=gff, nrBasePerSeg=nrbps)
    save(segScore, file=file.path(indir, sprintf("segScore-%d.rda", as.integer(nrbps))),
         compress=TRUE)
  }  
}
