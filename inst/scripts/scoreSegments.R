options(error=recover, warn=2)
library("tilingArray")
source("/homes/huber/madman/Rpacks/tilingArray/R/scoreSegments.R")

if(!exists("gff")) {
  cat("Loading probeanno.rda ")
  load("probeAnno.rda")
}

if(!exists("x")) {
  cat("x.rda\n")
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

chrs = 1:17
indirList = c("segmentation-3polyA", "seg-polyA-050428", "seg-tot-050421", "seg-polyA-050518")[4]
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
    ## segScore = addDirectHybe(segScore)
    cat("NOT ADDING DIRECT HYBE\n")
    save(segScore, file=file.path(indir, sprintf("segScore-%d.rda", as.integer(nrbps))),
         compress=TRUE)
  }  
}
