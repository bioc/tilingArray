## 5 plots:
## a) distribution of 5' UTR lengths
## b) distribution of 3' UTR lengths
## c) scatterplot 5' vs 3' UTR lengths
## d) scatterplot mean level vs 5' UTR length
## e) scatterplot mean level vs 3' UTR length
## The poly-A version is for the paper,
## the total RNA version is for the supplement.
##
## Criteria for good gene-containing segments:
## a. contains less than 50% (maxDuplicated) multiple-hit probes
## b. contains exactly one gene, fully (minOverlap): isOneGene 
## c. the gene is not spliced, and it is Validated or Uncharacterized (not Dubious): isGoodGene
## d. the mean level is above threshold
## e. on both flanks level is down
## f. moving average of 3 probes should not deviate too far to below

library("tilingArray")
library("prada")
source("colorRamp.R")
source("Figures/readSegments.R") 

graphics.off();
x11(width=14, height=9)
par(mfrow=c(2,3))
  
minOverlap=1
maxDuplicated=0.5

cols = brewer.pal(12, "Paired")

## For criterion c
goodGenes = gff$Name[gff$feature=="gene" &
  gff$orf_classification %in% c("Uncharacterized", "Verified")]

splicedGenes1 = sort(gff$Name[gff$feature=="intron"])
splicedGenes2 = names(which(table(gff$Name[gff$feature=="CDS"])>=2))
goodGenes = setdiff(goodGenes, union(splicedGenes1, splicedGenes2))

utr = vector(mode="list", length=length(rnaTypes))
names(utr)=rnaTypes
      
for(rt in rnaTypes) {
  s        = get("segScore", get(rt))
  isUnique = (s$frac.dup < maxDuplicated)        ## a
  isUnanno = (s$same.feature=="")
  thresh   = calcThreshold(s$level, sel=isUnique&isUnanno, main=rt)
  cat(rt, ": thresh=", signif(thresh, 2), "\n", sep="")

  isOneGene  = (listLen(strsplit(s$same.feature, split=", "))==1) & (s$same.overlap >= minOverlap) ## b
  isGoodGene = (s$same.feature %in% goodGenes)   ## c

  k = 2:(nrow(s)-1)  ## d+e
  isTranscribed = rep(FALSE, nrow(s))
  isTranscribed[k] = (s$level[k]>=thresh) & (s$level[k-1]<thresh) & (s$level[k+1]<thresh)

  sel = isUnique & isAnno & isGoodGene & isTranscribed

  utr5 = s$same.dist5[sel]
  utr3 = s$same.dist3[sel]

  z = cbind(utr5, utr3)
  rownames(z) = s$same.feature[sel]
  colnames(z) = c("5' UTR", "3' UTR")
  utr[[rt]] = z

  br=50
  hist(utr5[utr5<1000],
       col=cols[1], breaks=br,
       main=paste(rt, ": 5' UTR", " (", length(utr5), ")", sep=""))
  hist(utr3[utr3<1000], 
       col=cols[3], breaks=br, main=paste(rt, ": 3' UTR", sep=""))
  plot(utr5, utr3, pch=16, col=cols[2])
}

dev.copy(pdf, file="utrmap-hists.pdf", width=10, height=7); dev.off()

## common:
comm = intersect(rownames(utr[[1]]), rownames(utr[[2]]))
x11(width=9, height=5)
par(mfrow=c(1,2))
for(i in 1:2)
  smoothScatter(log(utr[[1]][comm,i]), log(utr[[2]][comm,i]),
       main=paste(colnames(utr[[1]])[i], " (", length(comm), ")", sep=""),
       xlab=rnaTypes[1], ylab=rnaTypes[2])

dev.copy(pdf, file="utrmap-scatter.pdf", width=8, height=4.8); dev.off()

out = file("utrmap.txt", open="wt")
cat("Selection criteria: segment must\n- contain exactly 1 annotated (",
    "'Verified' or 'Uncharacterized') gene.\n(There are ", length(goodGenes),
    " of these.)\n",
    "- less than ", signif(100*maxDuplicated, 2), "% of its sequence must ",
    "be hit by duplicated probes.\n",
    "- must be expressed above the 5% FDR threshold.\n\n", 
    sprintf("%10s  %6s  %5s\n", "RNA type", "side", "mode"), file=out, sep="")
for(rt in rnaTypes)
  for(i in 1:2)
    cat(sprintf("%10s  %6s  %5d\n", rt, colnames(utr[[rt]])[i],
                as.integer(shorth(utr[[rt]][comm, i]))), file=out)
close(out)
