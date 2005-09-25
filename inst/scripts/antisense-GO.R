##
## run this script after "tableSegments"
##
library("GOstats")
source(scriptsDir("GOHyperG.R"))

interact= TRUE

outfile = "antisense-GO"
if(!interact)
  sink(paste(outfile, "txt", sep="."))

if(!exists("utr"))
  load("utr.rda")

if(!exists("asCand")) {
  
  ## a matrix that corresponds to all ORFs
  nrGenes = length(featNames$"annotated ORFs")
  asCand = data.frame(
    "any overlap, filtered" = rep(FALSE, nrGenes),
    "3' overlap, filtered" = rep(FALSE, nrGenes),
    "5' overlap, filtered" = rep(FALSE, nrGenes),
    "any overlap, all" = rep(FALSE, nrGenes),
    "3' overlap, all" = rep(FALSE, nrGenes),
    "5' overlap, all" = rep(FALSE, nrGenes),
    "level"  = rep(as.numeric(NA), nrGenes),
    "3' UTR" = rep(as.numeric(NA), nrGenes),
    "5' UTR" = rep(as.numeric(NA), nrGenes),
    check.names=FALSE)
  
  rownames(asCand)=featNames$"annotated ORFs"
  
  ## UTR lengths
  whUTR = utr[["seg-polyA-050909"]]
  stopifnot(colnames(whUTR) %in% colnames(asCand))
  for(j in colnames(whUTR))
    asCand[rownames(whUTR), j] = whUTR[, j]

  ##
  ## Here we could use both poly-A and total RNA. In the paper we have poly-A
  ## only. The moral of the results doesn't change if we pool both
  ##
  allRnaTypes = rnaTypes[1]
  for(rt in allRnaTypes) {
    s = cs[[rt]]

    stopifnot("seg-polyA-050909" %in% allRnaTypes)
    if(rt == "seg-polyA-050909") {
      ## poly-A only....
      
      ## numeric vector, for each gene, the index of the segment it
      ##   is contained in (if any)
      inSegment = rep(as.numeric(NA), nrow(asCand))
      names(inSegment) = rownames(asCand)
      spFeatInSeg = strsplit(s$featureInSegment, ", ") 
      for(i in seq(along=spFeatInSeg)) {
        w = match(spFeatInSeg[[i]], names(inSegment))
        if(any(!is.na(inSegment[w])))
          cat("Warning - featureInSegment is defective - ", spFeatInSeg[[i]], "\n")
        inSegment[w] = i
      }
      asCand[, "level"] = s$level[inSegment]
    }
    
    
    ## numeric vector, for each gene, the index of the segments it
    ##   that list it as oppositeFeature
    oppositeSegment = vector(mode="list", length=nrow(asCand))
    names(oppositeSegment) = rownames(asCand)
    spOppoFeat  = strsplit(s$oppositeFeature, ", ")
    for(i in seq(along=spOppoFeat)) {
      w = match(spOppoFeat[[i]], names(oppositeSegment))
      for(v in w[!is.na(w)])
        oppositeSegment[[v]] = c(oppositeSegment[[v]], i)
    }
    
    cat("\n", rt, ":", sep="")
    for(j in 1:nrow(asCand)) {
      if(j%%1000==0)cat(j, "")
      
      ## segments that have this gene as their oppositeFeature
      iseg = oppositeSegment[[j]]
      if(length(iseg)>=1){
        
        wh1 = which(s$category[iseg] == "novel antisense - filtered")
        wh2 = c(wh1, which(s$category[iseg] == "novel antisense - unassigned"))
        
        if(length(wh2)>=1) {
          ## isGene is defined in tableSegments
          igff = which(gff$Name==rownames(asCand)[j] & isGene)
          
          asCand[j,  "any overlap, all"] = TRUE
          overl  = c(any( (s[iseg[wh2], "start"]-50) <= gff$start[igff] ),
            any( (s[iseg[wh2], "end"]  +50) >= gff$end[igff]   ))
          if(gff$strand[igff]=="-") ## reverse
            overl = rev(overl)
          asCand[j,  "3' overlap, all"] = overl[2]
          asCand[j,  "5' overlap, all"] = overl[1]
        }
      
        if(length(wh1)>=1) {
          asCand[j,  "any overlap, filtered"] = TRUE
          overl  = c(any( (s[iseg[wh1], "start"]-50) <= gff$start[igff] ),
            any( (s[iseg[wh1], "end"]  +50) >= gff$end[igff]   ))
          if(gff$strand[igff]=="-") ## reverse
            overl = rev(overl)
          asCand[j,  "3' overlap, filtered"] = overl[2]
          asCand[j,  "5' overlap, filtered"] = overl[1]
        }
        
      } ## if length(iseg)
    } ## for j
  } ## for it
} ## if(!exists("asCand"))

## 1. Is there a preference for overlap of antisense transcripts with 3' or 5' UTR?
cat("Is there a preference for overlap of antisense transcripts with 3' or 5' UTR?\n")
for(cn in colnames(asCand)[1:6])
  cat(sprintf("%30s: %4d\n", cn, as.integer(sum(asCand[, cn]))))
cat("\n\n")

## 2. Do genes with antisense transcripts have different levels?
g1 = which(asCand$"any overlap, filtered")
g2 = which(!asCand$"any overlap, all")
wt = wilcox.test(asCand$level[g1],  asCand$level[g2])
cat("Do genes with antisense transcripts have different levels from those without?\n")
cat(sprintf("Median of genes opposite 'novel antisense - filtered' segment: %8g\n",
            median(asCand$level[g1], na.rm=TRUE)),
    sprintf("Median of genes not opposite any 'novel antisense' segment: %8g\n",
            median(asCand$level[g2], na.rm=TRUE)),
    sprintf("Wilcoxon test p=%10s\n\n\n", format.pval(wt$p.value)), sep="")

## 3. Compare UTR lengths for UTRs that are overlapped by an antisense transcript
##  (filtered) versus those that are not (unfiltered/all)
cat("Compare UTR lengths for UTRs that are overlapped by an antisense transcript\n",
    "(filtered) versus those that are not (unfiltered/all):", sep="")
for(end in c("3'", "5'")) {
  whUTR = paste(end, "UTR")
  for(what in c("filtered", "all")) {
    g1 = which( asCand[, paste(end, "overlap,", what)] & !is.na(asCand[, whUTR]))
    g2 = which(!asCand[, "any overlap, all"] & !is.na(asCand[, whUTR]))
    x1 = asCand[g1, whUTR]
    x2 = asCand[g2, whUTR]
    wt = wilcox.test(x1,  x2)
    cat(sprintf("\n%2s, using >>%s<< putative antisense segments:\n", end, what),
        sprintf("Median length of UTRs of %d genes opposite a 'novel antisense' segment: %8g\n",
                length(x1), median(x1)),
        sprintf("Median length for %d genes not opposite any 'novel antisense' segment, and with UTR length data: %8g\n",
                length(x2), median(x2)),
        sprintf("Wilcoxon test p=%10s\n", format.pval(wt$p.value)), sep="")
  }
}


stop()

#######################################################################
##  look for enriched GO classes - use combination of poly-A and total
#######################################################################

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


what = "all"
cat("\n\n\nGO category analysis (", what, "):\n", sep="")
cat("======================================\n\n")
GOHyperG(rownames(asMat)[asMat[, what]])

if(!interact) {
  dev.copy(pdf, paste(outfile, what, "pdf", sep="."), width=12, height=6.3)
  dev.off()
}

if(!interact)
  sink()

