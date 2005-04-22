library("tilingArray")
library("prada"); source("colorRamp.R")
source("Figures/readSegments.R") 

graphics.off();
x11(width=14, height=9)
par(mfrow=c(2,3))
  
minOverlap=1
maxDuplicated=0.5

cols = brewer.pal(12, "Paired")

goodGenes = gff$Name[gff$feature=="gene" &
  gff$orf_classification %in% c("Uncharacterized", "Verified")]

utr = vector(mode="list", length=length(rnaTypes))
names(utr)=rnaTypes
      
for(rt in rnaTypes) {
  s        = get("segScore", get(rt))
  isUnique = (s$frac.dup < maxDuplicated) 
  isUnanno = (s$same.feature=="")
  thresh   = calcThreshold(s$level, sel=isUnique&isUnanno, main=rt)
  cat(rt, ": thresh=", signif(thresh, 2), "\n", sep="")

  isAnno     = (listLen(strsplit(s$same.feature, split=", "))==1) & (s$same.overlap >= minOverlap)
  isGoodGene = s$same.feature %in% goodGenes
  isTranscribed = (s$level>=thresh)

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
