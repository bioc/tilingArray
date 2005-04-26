options(error=recover, warn=2)
library("tilingArray")
## source("/homes/huber/madman/Rpacks/tilingArray/R/scoreSegments.R")

if(!exists("gff"))
  load("probeAnno.rda")

chrs = 1:16

for(indir in c("segmentation-3polyA", "seg-tot-050421")) {
  cat(indir, "\n")

  if(TRUE) {
    s  = new.env()
    cat("Loading ")
    for(chr in chrs) {
      for(strand in c("+", "-")) {
        fn = file.path(indir, paste(chr, strand, "rda", sep="."))
        cat(chr, ".", strand, " ", sep="")
        load(fn)
        assign(paste(chr, strand, "seg", sep="."), seg, envir=s)
        assign(paste(chr, strand, "dat", sep="."), dat, envir=s)
      }
    } ## for chr
    cat("\n")
  } else {
    cat("NOT LOADING DATA FILES!\n")
  }
  
  ## For the definition of pseudogenes at SGD, see Docs/PseudogenesAtSGD.pdf
  for(nrbps in c(1000, 1500, 2000)) {
    cat(">> ", nrbps, "<<\n")
    segScore = scoreSegments(s, gff=gff, nrBasePerSeg=nrbps)
    save(segScore, file=file.path(indir, sprintf("segScore-%d.rda", as.integer(nrbps))),
         compress=TRUE)
  }  
}
