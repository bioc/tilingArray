options(error=recover, warn=2)
library("tilingArray")

if(!exists("gff"))
  load("probeAnno.rda")

indir = "segmentation-050209v4"
chrs = 1:17

if(!exists("s")) {
  s  = new.env()
  totcp = 0
  cat("Loading ")
  for(chr in chrs) {
    for(strand in c("+", "-")) {
      fn = file.path(indir, paste(chr, strand, "rda", sep="."))
      cat(chr, ".", strand, " ", sep="")
      load(fn)
      assign(paste(chr, strand, "seg", sep="."), seg, envir=s)
      assign(paste(chr, strand, "dat", sep="."), dat, envir=s)
      totcp = totcp + round(max(dat$x)/nrBasePerSeg)
    }
  } ## for chr
  cat("\n")
} ## if

segScore = scoreSegments(s, gff=gff)

save(segScore, file=file.path(indir, "segScore.rda"), compress=TRUE)
