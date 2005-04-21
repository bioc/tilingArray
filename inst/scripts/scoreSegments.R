options(error=recover, warn=2)
library("tilingArray")
source("/home/huber/madman/Rpacks/tilingArray/R/scoreSegments.R")

indir = c("segmentation-3polyA", "seg-tot-050418")[2]
cat(indir, "\n")

if(!exists("gff"))
  load("probeAnno.rda")

chrs = 1:17
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

## For the definition of pseudogenes at SGD, see Docs/PseudogenesAtSGD.pdf
## for(nrbps in c(1500, 1750, 2000, 2250)) {
for(nrbps in c(2000)) {
  segScore = scoreSegments(s, gff=gff, nrBasePerSeg=nrbps)
  save(segScore, file=file.path(indir, sprintf("segScore-%d.rda", as.integer(nrbps))),
                   compress=TRUE)
}  
