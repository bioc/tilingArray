library("tilingArray")
source("scripts/readSegments.R") 
source("scripts/calcThreshold.R") 
source("scripts/categorizeSegments.R") 
source("scripts/writeSegmentTable.R")

options(error=recover, warn=0)

interact = F
what=c("pie", "wst", "length", "lvsx", "cons", "conswex")

if(!interact & exists("tab"))
  rm(tab)

n = length(rnaTypes)

if(!interact)
  sink("tableSegments.txt")

if(!exists("tab")) {
  graphics.off()
  if(interact)
    x11(width=10, height=n*3)
  else
    pdf(width=11, height=n*4)
  par(mfrow=c(n,1))
  
  tab = vector(mode="list", length=length(rnaTypes))
  names(tab)=rnaTypes
  
  cat("Categorizing segments:\n",
      "======================\n", sep="")
  for(rt in rnaTypes) {
    s = get("segScore", get(rt))
    tab[[rt]] = categorizeSegmentsPie(s)
  }

  if(!interact)
    dev.off()
}


colors = c(brewer.pal(9, "Pastel1")[c(2:6, 1)])
names(colors) =c("verified", "ncRNA", "uncharacterized", "dubious",
     "unA", "unI")

##
## PIE
##
if("pie" %in% what){
  if(!interact)
    pdf("tableSegments-pie.pdf", width=7*n, height=4.8)
  par(mfrow=c(1,n))
  cat("Counts for the pie chart:\n",
      "=========================\n",
      "For known features, counts are unique IDs. While in typical cases\n",
      "there is a 1:1 mapping between segments and IDs, this is not always\n",
      "the case: there can be several IDs per segment, and several segments\n",
      "per ID.\n",
      "For new segments, counts are non-consecutive blocks of one or several\n",
      "expressed segments flanked by a non-expressed segment on each side.\n",
      sep="")
  for(rt in rnaTypes) {
    ct  = tab[[rt]]$count
    
    ctp  = rbind(ct,  "total" = colSums(ct, na.rm=TRUE))
    perct = 100*ctp[,1]/ctp[,2]
    ctp  = cbind(ctp, "percent" = round(perct,1))
    rownames(ctp) = sub("unA", "unannot. (pot. antisense)", rownames(ctp))
    rownames(ctp) = sub("unI", "unannot. (isolated)", rownames(ctp))
  
    cat("\n", rt, "\n-----\n", sep="")
    print(ctp)
    
    px = ct[, "observed"]
    names(px) = rownames(ct)
    px = px[c("verified",  "ncRNA", "uncharacterized", "dubious", "unA", "unI")]
    pie(px, radius=0.75, main=longNames[rt],
        col=colors,
        labels=paste(names(px), " (", px, ")", sep=""))
  }
  if(!interact)
    dev.off()
}

##
## WRITE THE SEGMENT TABLE
##
if("wst" %in% what){
  for(rt in rnaTypes) {
    s = get("segScore", get(rt))
    sel = tab[[rt]]$category %in% c("verified", "uncharacterized", "ncRNA", "unA", "unI")
    s = cbind(category=tab[[rt]]$category, s)[sel, ]
    fn = file.path(indir[rt], "viz", "index.html")
    cat("Writing", fn, "\n")
    writeSegmentTable(s, title=longNames[rt], fn=fn)
  }
}

##
## LENGTH DISTRIBUTIONS
##
maxlen=5000
if("length" %in% what){
  if(!interact)
    pdf("tableSegments-lengths.pdf", width=14, height=n*3)
  par(mfrow=c(n,5))
  br = seq(0, maxlen, by=200)
  for(rt in rnaTypes) {
    s  = get("segScore", get(rt))
    ct = tab[[rt]]$category
    for(lev in c("verified",  "uncharacterized", "ncRNA", "unA", "unI")) {
      len = s$length[ct == lev]
      len[len>maxlen]=maxlen
      hist(len, breaks=br, col=colors[lev], main=paste(rt, lev))
    }
  }
  if(!interact)
    dev.off()
}

##
## LENGTH VERSUS EXPRESSION LEVEL
##
if("lvsx" %in% what){
  if(!interact) {
    pdf("tableSegments-lvsx.pdf", width=14, height=n*3)
    pch="."
  } else {
    pch=18
  }
  par(mfrow=c(n, 2))
  maxlen=5000
  br = seq(0, maxlen, by=200)
  for(rt in rnaTypes) {
    s  = get("segScore", get(rt))
    ct = tab[[rt]]$category
    levs = c("verified", "unI")
    ylim = quantile(s$level[ct %in% levs], probs=c(0.01, 0.99), na.rm=TRUE)
    for(lev in levs) {
      len = s$length[ct == lev]
      len[len>maxlen]=maxlen
      exl = s$level[ct == lev]
      plot(len, exl, pch=pch, ylim=ylim, main=paste(longNames[rt], ": ", lev, sep=""),
           ylab="expression level", xlab="length")
      lf = loess(exl ~ len)
      slen = sort(len)
      lines(slen, predict(lf, newdata=slen), col="blue")
      ## smoothScatter(len, exl, ylim=ylim, pch=pch, main=paste(rt, lev))
    }
  }
  if(!interact)
    dev.off()
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
    
    alignableLen = numeric(nrow(s))
    alignableLen[  br[[1]][theBest] ] = numNucMatch[theBest]
    
    hit[,b] = sapply(sp, function(segments) {
      sum(alignableLen[segments]) / sum(s$length[segments]) 
    })
  }
  data.frame(number=listLen(sp), fraction=rowMeans(hit))
}

##
## cons
##

if("cons" %in% what){
  cat("\n\nFraction of alignable sequence (percent):\n",
      "=========================================\n",
      "Here, 'number' is the number of segments.\n", sep="")
  
  for(rt in rnaTypes) {
    s  = get("segScore", get(rt))
    ct = tab[[rt]]$category
    sp = split(seq(along=ct), ct)
    selected = c("verified", "uncharacterized" , "ncRNA",  "unA", "unI")
    stopifnot(all(selected %in% names(sp)))
    sp = sp[selected]
    sp$unexpressed=which(s$featureInSegment=="" & s$segmentInFeature=="" &
        s$oppositeFeature=="" & s$level <= quantile(s$level, 0.2, na.rm=TRUE))
    sp$"whole genome" = seq(along=ct)
    
    hit = calchit(sp, rt)
    
    cat("\n", rt, "\n-----\n", sep="")
    print(round(hit,1))
  }
}

##
## conswex
##
if("conswex" %in% what){
  nrlevs = 5
  minlen    = 200
  cat("\n\nGrouping of conservation scores by length and expression levels\n",
      "using ", nrlevs,
      " quantile groups of equal size, respectively.\n",
      "For 'expression', scores are calculated for segments with length >= ", minlen, ".\n",
      "==========================================================================\n", 
      sep="")
  if(!interact)
    pdf("tableSegments-conswex.pdf", width=4*n, height=8)
  par(mfcol=c(2,n))
  cols = c("#303030", "#0000e0")
  pchs = c(15,16)
  
  segClasses = c("verified", "unI")
  segClassesLong = c("verified genes", "unannotated, isolated")
  fraction = matrix(NA, nrow=nrlevs, ncol=length(segClasses))
  colnames(fraction) = segClasses
  
  for(rt in rnaTypes) {
    s   = get("segScore", get(rt))
    ct  = tab[[rt]]$category
    cat("\n", rt, "\n-----\n", sep="")

    for(ev in c("expression", "length")) {
      for(segClass in segClasses) {
        wh = which(ct == segClass)
        v  = switch(ev,
          "expression" = s$level[wh],
          "length" =     s$length[wh],
          stop("Zapperlot"))
        br = quantile(v, (0:nrlevs)/nrlevs)
        sp = split(wh, cut(v, breaks=br))
        names(sp) = switch(ev,
          "expression" = paste("<=", round(br[-1], 1)),
          "length"     = paste("<=", as.integer(br[-1], 0)))
        
        if( ev=="expression")
          sp = lapply(sp, function(i) i[s$length[i]>=minlen])

        hit = calchit(sp, rt)
        fraction[, segClass] = hit$fraction
        cat("\nby '", ev, "', for '", segClass, "' segments\n", sep="")
        print(round(hit,1))
      }
      matplot(fraction, xaxt="n", type="b", main=paste(ev, " (", rt, ")", sep=""),
              lty=1, lwd=2, pch=pchs, col=cols,
              ylab="Fraction of alignable sequence", xlab=ev,
              ylim=c(0,100))
      axis(side=1, at = 1:nrlevs, labels = names(sp))
      if(rt=="polyA2" & ev=="expression")
        legend(x=1, y=100, legend=segClassesLong, lty=1, lwd=2, pch=pchs, col=cols)
    }
  }
  if(!interact)
    dev.off()
}

if(!interact)
  sink()
