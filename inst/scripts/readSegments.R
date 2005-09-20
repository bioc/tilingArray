if(!exists("gff")) {
  cat("Loading probeAnno.rda\n")
  load("probeAnno.rda")
}

longNames = c("seg-polyA-050811"    = "poly-A RNA",
              "seg-tot-050811"      = "total RNA",
              "seg-dir-050811"      = "direct",
              "seg-odT-050811"      = "oligo-dT",
              "seg-polyA0420-050811"= "clean poly-A",
              "seg-polyA-050909"    = "poly-A RNA",
              "seg-tot-050909"      = "total RNA",
              "seg-dir-050909"      = "direct",
              "seg-odT-050909"      = "oligo-dT",
              "seg-polyA0420-050909"= "clean poly-A")


rtOK = rnaTypes %in% names(longNames)
if(!all(rtOK))
  stop(paste("'longNames' not defined for: '",
       paste(rnaTypes[!rtOK], collapse="', '"),
       "'.", sep=""))

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


