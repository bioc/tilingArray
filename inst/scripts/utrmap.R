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
library("geneplotter")

interact=(!TRUE)
options(error=recover, warn=0)
graphics.off()

## source("colorRamp.R")
source("scripts/readSegments.R")
source("scripts/calcThreshold.R") 
source("scripts/categorizeSegments.R") 
source("scripts/writeSegmentTable.R")


fdim = c(10.5, 7)
if(interact) {
  x11(width = fdim[1], height = fdim[2])
  pch=16
} else {
  sink("utrmap.txt")
  pch="."
}
  
cols = brewer.pal(12, "Paired")
trsf = function(x) log(x+1, 10)
## trsf = function(x) sqrt(x)

investigateExpressionVersusLength = function(lev, len, main) {
  theCut = cut(len, breaks=quantile(len, probs=c(0, 0.8, 0.95, 1)))
  ## e = tapply(len, cut(lev, breaks=3), ecdf)
  e = tapply(lev, theCut, ecdf)
  theCol = cols[1]
  plot(e[[1]], xlab="level", main=main)
  for(i in 2:length(e)) {
    theCol = cols[i*2]
    lines(e[[i]], col.hor=theCol, col.points=theCol, col.vert=theCol)
  }
}

utr = vector(mode="list", length=length(rnaTypes))
names(utr)=rnaTypes
      
for(rt in rnaTypes) {
  cat("\n", rt, "\n-----\n", sep="")
  s = categorizeSegmentsUTRmap(get(rt))
  s = s[!is.na(s$goodUTR), ]
  
  ##
  ## WRITE THE SEGMENT TABLE
  ##
  fn = file.path(indir[rt], "viz", "utrmap.html")
  if(TRUE){
    cat("Writing", nrow(s), "UTRs to", fn, "\n")
    writeSegmentTable(s, title=paste(nrow(s), "UTR maps from", longNames[rt]), fn=fn,
                      sortBy = "goodUTR", sortDecreasing=TRUE)
  } else {
    cat(">>> Attention: not writing", fn, "<<<\n") 
  }

  cat("5' UTR length distribution summary:\n")
  print(summary(s$utr5))
  cat("3' UTR length distribution summary:\n")
  print(summary(s$utr3))
  
  z = cbind(s$utr5, s$utr3)
  rownames(z) = s$geneInSegment
  colnames(z) = c("5' UTR", "3' UTR")
  utr[[rt]] = z

  ##
  ## LENGTH HISTOGRAM
  ##
  if(!interact) {
    pdf(file=paste("utrmap-", rt, ".pdf", sep=""), width = fdim[1], height = fdim[2])
  }
  par(mfrow = c(2, 3))
  br=50
  hist(s$utr5[s$utr5<1000],
       col=cols[1], breaks=br, xlab="length of 5' UTR",
       main=paste(longNames[rt], ": 5' UTR", " (", length(s$utr5), ")", sep=""))
  hist(s$utr3[s$utr3<1000], xlab="length of 3' UTR",
       col=cols[3], breaks=br, main=paste(longNames[rt], ": 3' UTR", sep=""))

  ##
  ## GET THE CDS LENGTHS AS WELL
  ##
  mt = match(s[,"geneInSegment"], gff[,"Name"])
  stopifnot(!any(is.na(mt)))
  cdslen = gff[mt, "end"]-gff[mt, "start"]

  investigateExpressionVersusLength(s$level, s$utr3, "length of 3' UTR")
  investigateExpressionVersusLength(s$level, s$utr5, "length of 5' UTR")
  investigateExpressionVersusLength(s$level, cdslen, "length of CDS")

  if(!interact)
    dev.off()
}
cat("\n\n")


## common:
comUTR = intersect(rownames(utr[[1]]), rownames(utr[[2]]))
## union:
allUTR = union(rownames(utr[[1]]), rownames(utr[[2]]))
for(j in seq(along=utr))
  cat(nrow(utr[[j]]), " UTRs for ", rnaTypes[j], ", ", sep="")
cat("\n", length(comUTR)," in both ", paste(rnaTypes, collapse=" and "), ", ", 
    length(allUTR), " altogether.\n", sep="")

if(interact) {
  x11(width=8, height=4)
} else {
  pdf(file=paste("utrmap-scatter.pdf", sep=""), width=8, height=4)
}

par(mfrow=c(1,2))
for(i in 1:2){
  smoothScatter(trsf(utr[[1]][comUTR,i]), trsf(utr[[2]][comUTR,i]),
       main=paste("log10 of length of ", colnames(utr[[1]])[i], " (", length(comUTR), ")", sep=""),
       xlab=longNames[rnaTypes[1]], ylab=longNames[rnaTypes[2]])
  abline(a=0, b=1, col="#606060")

  
}

w = (abs(utr[[1]][comUTR, 1] - utr[[2]][comUTR, 1]) < 10 &
     abs(utr[[1]][comUTR, 2] - utr[[2]][comUTR, 2]) < 10 )
u5 = utr[[1]][comUTR[w], 1]
u3 = utr[[1]][comUTR[w], 2]

if(!interact) {
  sink()
  dev.off()
}

