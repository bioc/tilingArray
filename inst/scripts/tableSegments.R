library("tilingArray")
source("scripts/readSegments.R") 
source("scripts/calcThreshold.R") 
source("scripts/categorizeSegments.R") 

options(error=recover, warn=0)

doSink=TRUE
what=c("pie", "length", "cons", "conswex")[1:4]

if(doSink)
  sink("tableSegments.txt")

if(!exists("tab")) {
  graphics.off(); x11(width=9, height=5)
  par(mfrow=c(2,1))
  
  tab = vector(mode="list", length=length(rnaTypes))
  names(tab)=rnaTypes
  
  for(rt in rnaTypes) {
    s = get("segScore", get(rt))
    tab[[rt]] = categorizeSegmentsPie(s)
  }

  dev.copy(pdf, "tableSegments-thresh.pdf", width=11, height=8); dev.off()
}


colors = c(brewer.pal(9, "Pastel1")[c(2:6, 1)])
names(colors) =c("verified", "ncRNA", "uncharacterized", "dubious",
     "unA", "unI")

##
## PIE
##
if("pie" %in% what){
  par(mfrow=c(1,2))
  cat("Counts for the pie chart:\n",
      "=========================\n",
      "For known features, counts are unique IDs.\n",
      "For new segments, counts are non-consecutive blocks of 1..n segments\n",
      sep="")
  for(rt in rnaTypes) {
    ct  = tab[[rt]]$count
    
    perct = 100*ct[,1]/ct[,2]
    ctp  = cbind(ct, "percent" = round(perct,1))
    ctp  = rbind(ctp,  "total" = colSums(ctp, na.rm=TRUE))
    
    cat("\n", rt, "\n-----\n", sep="")
    print(ctp)
    
    px = ct[, "observed"]
    names(px) = rownames(ct)
    px = px[c(1,3,2,4,6,5)]
    pie(px, radius=0.75, main=c(polyA="poly-A RNA", tot="total RNA")[rt],
        col=colors,
        labels=paste(names(px), " (", px, ")", sep=""))
  }
  dev.copy(pdf, "tableSegments-pie.pdf", width=14, height=4.8); dev.off()
}

##
## LENGTH DISTRIBUTIONS
##
if("length" %in% what){
  par(mfrow=c(2,5))
  maxlen=5000
  br = seq(0, maxlen, by=200)
  for(rt in rnaTypes) {
    s  = get("segScore", get(rt))
    ct = tab[[rt]]$category
    for(lev in levels(ct)[c(1:2, 4:6)]) {
      len = s$length[ct == lev]
      len[len>maxlen]=maxlen
      hist(len, breaks=br, col=colors[lev], main=paste(rt, lev))
    }
  }
  dev.copy(pdf, "tableSegments-lengths.pdf", width=14, height=6); dev.off()
}

##
## CONSERVATION
##

if(!exists("blastres")) {
  blastResultFiles = c("Sbay_contigs.out", "Smik_contigs.out",
    "Spar_contigs.out")
  names(blastResultFiles) = c("S.bayanus", "S.mikatae",
         "S.paradoxus")
  names(rnaTypes)=rnaTypes
  blastres = lapply(rnaTypes, function(rt)
    lapply(blastResultFiles, function(f)
	 read.table(file.path(indir[rt], "fasta", f),
              sep="\t", as.is=TRUE, header=FALSE)))
}
    
calchit = function(sp, rt) {  
  hit = matrix(NA, nrow=length(sp), ncol=length(blastResultFiles))
  rownames(hit) = names(sp)
  colnames(hit) = names(blastResultFiles)
  for(b in 1:ncol(hit)) {
    br = blastres[[rt]][[b]]
    
    ## numNucMatch = br[[3]]*br[[4]] ## 3=Percent identity, 4=Alignment length
    numNucMatch = 100*br[[4]] ## 3=Percent identity, 4=Alignment length

    ## split by name of query sequence (1) and just keep the hit with the
    ## highest value of numNucMatch
    spbyq = split(1:nrow(br), br[[1]]) ## 1=Identity of query sequence
    theBest = sapply(spbyq, function(i) i[which.max(numNucMatch[i])])
    
    prcid = numeric(nrow(s))
    prcid[  br[[1]][theBest] ] = numNucMatch[theBest]
    prcid = prcid / s$length
    
    hit[,b] = sapply(sp, function(segments) {
      mean( prcid[segments] )
    })
  }
  hit = data.frame(number=listLen(sp), fraction=rowMeans(hit))
  cat("\n", rt, "\n-----\n", sep="")
  print(round(hit,1))
  hit  
}


if("cons" %in% what){
  cat("\n\nFraction of alignable sequence (percent):\n",
      "=========================================\n",
      sep="")
  
  for(rt in rnaTypes) {
    s  = get("segScore", get(rt))
    ct = tab[[rt]]$category
    sp = split(seq(along=ct), ct)

    stopifnot(length(sp)==7)
    sp[[7]] = which(s$same.feature=="" & s$same.dist5 > 100 & s$same.dist3 > 100 &
                 s$level <= quantile(s$level, 0.2, na.rm=TRUE))
    sp[[8]] = seq(along=ct)
    names(sp)[8] = "whole genome"
    sp = lapply(sp, function(i) i[s$length[i]>=800 & s$length[i]<=1200])
  
    hit = calchit(sp, rt)
  }
 }


## does the hit rate depend on expression level?
if("conswex" %in% what){
  cat("\n\nanI, grouped by expression levels\n(in 5 quantile groups of equal size):\n",
      "=====================================\n", 
      sep="")
  par(mfrow=c(1,2))
  cols = 1:4

  segClasses = c("verified", "unI")
  segClassesLong = c("verified genes", "unannotated, isolated")
  nrexplevs = 5
  fraction = matrix(NA, nrow=nrexplevs, ncol=length(segClasses))
  colnames(fraction) = segClasses
  
  for(rt in rnaTypes) {
    s   = get("segScore", get(rt))
    ct  = tab[[rt]]$category
    for(segClass in segClasses) {
      wh = which(ct == segClass)
      sp = split(wh, cut(s$level[wh], breaks=quantile(s$level[wh], (0:nrexplevs)/nrexplevs)))
      ## sp = lapply(sp, function(i) i[s$length[i]>=700 & s$length[i]<=1400])
      fraction[, segClass] = calchit(sp, rt)$fraction
    }
    matplot(fraction, type="b", main=rt, lty=1, 
            ylab="Fraction of alignable sequence", xlab="expression",
            ylim=c(0,100), pch=16)
    if(rt=="polyA")
      legend(x=1, y=100, legend=colnames(hit)[-1], lty=1, pch=16, col=cols)
  }
  dev.copy(pdf, "tableSegments-conswex.pdf", width=12, height=6); dev.off()
}

if(doSink)
  sink()
