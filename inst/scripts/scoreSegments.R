options(error=recover, warn=2)
library("tilingArray")

if(!exists("gff"))
  load("probeAnno.rda")

indir = "segmentation-050209v4"
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
segScore = scoreSegments(s, gff=gff)

save(segScore, file=file.path(indir, "segScore.rda"), compress=TRUE)
    
