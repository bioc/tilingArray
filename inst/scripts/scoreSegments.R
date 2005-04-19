options(error=recover, warn=2)
library("tilingArray")
source("/home/huber/madman/Rpacks/tilingArray/R/scoreSegments.R")

if(!exists("gff"))
  load("probeAnno.rda")

indir = "segmentation-3polyA"
chrs = 1:17

if(!exists("s")) {
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
} ## if

## For the definition of pseudogenes at SGD, see Docs/PseudogenesAtSGD.pdf

for(nrbps in c(1500, 1750, 2000, 2250)) {
  segScore = scoreSegments(s, gff=gff, nrBasePerSeg=nrbps)
  save(segScore, file=file.path(indir, sprintf("segScore-%d-NEW.rda", as.integer(nrbps))),
                   compress=TRUE)
}  
