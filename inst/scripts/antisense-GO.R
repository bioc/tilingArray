##
## run this script after "tableSegments"
##
##  Associations between presence of antisense transcripts and GO categories
##
library("GOstats")
source(scriptsDir("GOHyperG.R"))

interact = !TRUE

outfile = "antisense-GO"
if(!interact)
  sink(paste(outfile, "txt", sep="."))

##
## use combination of poly-A and total
##
asMat = matrix(FALSE, nrow=length(featNames$"annotated ORFs"), ncol=2)
rownames(asMat)=featNames$"annotated ORFs"
colnames(asMat)=c("filtered", "all")

for(what in colnames(asMat)) {
  catgSel = list(filtered = "novel antisense - filtered",
                 all = c("novel antisense - filtered", "novel antisense - unassigned"))[[what]]
  for(rt in rnaTypes) {
    s = cs[[rt]]
    selSeg = which(s[, "category"] %in% catgSel)
    asGenes = unique(unlist(strsplit(s[selSeg, "oppositeFeature"], split=", ")))
    asGenes = intersect(asGenes, rownames(asMat))
    asMat[ asGenes, what ] = TRUE
  } ## for rt
} ## for what


if(!interact)
  pdf(paste(outfile, "pdf", sep="."), width=12, height=12)

par(mfrow=c(2,2))

for(what in colnames(asMat)) {
  cands = rownames(asMat)[asMat[, what]]
  cat("\nGO category analysis (", length(cands), " genes, from: '", what, "'):\n", sep="")
  cat("===================================================\n\n")
  GOHyperG(cands, plottitle=what,
           outtable=paste(outfile, "-table-", what, ".txt", sep=""))
}

if(!interact)
  dev.off()

if(!interact)
  sink()

