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
## d. the mean level is above threshold
## e. on both flanks level is down
## f. moving average of 3 probes should not deviate too far to below

library("tilingArray")
## source("colorRamp.R")
source("scripts/readSegments.R") 
source("scripts/categorizeSegments.R") 
source("scripts/writeSegmentTable.R")

graphics.off();
x11(width=14, height=9)
par(mfrow=c(2,3))
  
cols = brewer.pal(12, "Paired")

utr = vector(mode="list", length=length(rnaTypes))
names(utr)=rnaTypes
      
for(rt in rnaTypes) {
  s = get("segScore", get(rt))
  s$category = categorizeSegmentsUTRmap(s)
 
  ##
  ## WRITE THE SEGMENT TABLE
  ##
  s$gene = character(nrow(s))
  selgff = which(gff$feature=="gene")
  mt = match(s$featureInSegment, gff$Name[selgff])
  hasMatch    = !is.na(mt)
  hasGeneName = !is.na(gff$gene[selgff][mt])
  sel = hasMatch&hasGeneName
  s$gene[sel] = gff$gene[selgff][mt[sel]]
  
  wh  = which(!is.na(s$category))
  ord = order(s$category[wh])
  s = s[wh[ord], c( "category", "excurse", "sdLeft", "sdThis", "sdRight",
    "gene", "chr", "strand", "start", "end", "length", "level", "utr5", "utr3",
    "featureInSegment", "segmentInFeature", "oppositeFeature",  "frac.dup")]
  colnames(s)[1] = "score"
  s = cbind(rank=1:nrow(s), s)
  
  fn = file.path(indir[rt], "viz", "utrmap.html")
  cat("Writing", fn, "\n")

  writeSegmentTable(s, title=paste(nrow(s), "UTR maps from", longNames[rt]), fn=fn)

  browser()
  utr5 = s$same.dist5[ wh ]
  utr3 = s$same.dist3[ wh ]

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
