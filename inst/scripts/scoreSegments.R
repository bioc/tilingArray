### 2005-08-26: run function 'scoreSegments' on data; see tilingArray/R/scoreSegments.R for the used functions.

library("tilingArray")

# where to find new versions of the functions without reinstalling tilingArray:
rfuncDir <- "/ebi/research/huber/users/joern/tilingArray/R"
# where to find the scripts:
scriptsDir <-  "/ebi/research/huber/users/joern/tilingArray/inst/scripts"

source(file.path(rfuncDir, "scoreSegments.R"))

rnaTypes  = c("seg-polyA-050804",
  "seg-polyA-050810", "seg-polyA-050811", "seg-tot-050811",
  "seg-dir-050811" , "seg-odT-050811", "seg-polyA0420-050811")[-c(1,2)]

doNotLoadSegScore=TRUE
source(file.path(scriptsDir, "readSegments.R"))

names(indir) = indir = rnaTypes

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
    sel = ind[(uni==0) & (sta >= s[i, "start"]) & (end <= s[i, "end"])]
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

gffBS <- grep("binding_site", gff$feature)
table(gff[gffBS,"strand"])
#   -    .    +
#  15    0 3369

# compare with:
table(gff[grep("CDS$", gff$feature),"strand"])
#    -    .    +
# 3049    0 3177

