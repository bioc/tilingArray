if(!exists("gff")) {
  cat("Loading probeAnno.rda\n")
  load("probeAnno.rda")
}

longNames = c("seg-polyA-050525"="poly-A RNA",
  "seg-tot-050525"="total RNA",
  "seg-tot2-050525"="total RNA (v2)")

names(indir) = indir = rnaTypes

for(rt in rnaTypes) {
  if(!exists(rt)) {
    assign(rt, new.env())
    chrs = 1:17
    cat("Loading", rt, ": ")
    for(chr in chrs) {
      cat(chr, "")
      for(strand in c("+", "-")) {
        fn = paste(chr, strand, "rda", sep=".")
        load(file.path(indir[rt], fn))  
        assign(paste(chr, strand, "seg", sep="."), seg, envir=get(rt))
        assign(paste(chr, strand, "dat", sep="."), dat, envir=get(rt))
      }
    } ## for chr
    if(!exists("doNotLoadSegScore")) {
      fn="segScore-1500.rda"
      cat(fn)
      load(file.path(indir[rt], fn), envir=get(rt))
    }
    cat("\n")
  } ## if
} ## rnaTypes


