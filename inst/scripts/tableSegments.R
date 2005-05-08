library("tilingArray")

graphics.off()
options(error=recover, warn=2)
interact = FALSE
what     = c("pie", "wst", "length", "lvsx", "cons")

if(!interact)
  sink("tableSegments.txt")

source("scripts/readSegments.R") 
source("scripts/calcThreshold.R") 
source("scripts/categorizeSegments.R") 
source("scripts/writeSegmentTable.R")

##
## CATEGORIZE
##
cs = vector(mode="list", length=length(rnaTypes))
names(cs)=rnaTypes

cat("Categorizing segments:\n",
    "======================\n", sep="")
for(rt in rnaTypes) {
  cs[[rt]] = categorizeSegmentsPie(get(rt))
}



fillColors = c(brewer.pal(9, "Pastel1")[c(2:6, 1)])
names(fillColors) =c("verified", "ncRNA", "uncharacterized", "dubious",
      "unAnti", "unIso")

lineColors = c(brewer.pal(9, "Set1"))[c(1, 5, 4, 3, 2, 9, 7)]
names(lineColors) =c("verified", "uncharacterized", "ncRNA", 
                 "unAnti", "unIso", "unannotated, unexpressed", "whole genome")

##
## PIE
##
if("pie" %in% what){
  if(!interact)
    pdf("tableSegments-pie.pdf", width=7*length(rnaTypes), height=4.8)
  par(mfrow=c(1, length(rnaTypes)))
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
    cat("\n", rt, "\n-----\n", sep="")
    selectedCategories = c("verified", "ncRNA", "uncharacterized", "dubious", "unIso", "unAnti")
    stopifnot(all(selectedCategories %in% names(fillColors)))
    ct  = cs[[rt]]$count
    stopifnot(all(selectedCategories %in% rownames(ct)))

    print(cbind(ct, "percent" = round(100*ct[,1]/ct[,2], 1)))

    px = ct[selectedCategories, "observed"]
    pie(px, radius=0.75, main=longNames[rt],
        col=fillColors[selectedCategories],
        labels=paste(selectedCategories, " (", px, ")", sep=""))
    rm(list=c("ct", "px", "selectedCategories"))
  }
  if(!interact)
    dev.off()
}

##
## WRITE THE SEGMENT TABLE
##
if("wst" %in% what){
  cat("\n\n")
  selectedCategories = c("verified", "uncharacterized", "ncRNA", "unAnti", "unIso", "unDubious")
  for(rt in rnaTypes) {
    s   = cs[[rt]]$s
    stopifnot(all(selectedCategories %in% levels(s$category)))
    sel = (s$category %in% selectedCategories)
    fn  = file.path(indir[rt], "viz", "index.html")
    writeSegmentTable(s[sel, ], title=longNames[rt], fn=fn)
    rm(list=c("s", "sel", "fn"))
  }
  cat("\n\n")
}

##
## LENGTH DISTRIBUTIONS
##
maxlen=5000
if("length" %in% what){
  selectedCategories = c("verified", "uncharacterized", "ncRNA", "unAnti", "unIso")
  stopifnot(all(selectedCategories %in% names(fillColors)))
  if(!interact)
    pdf("tableSegments-lengths.pdf", width=14, height=length(rnaTypes)*3)
  par(mfrow=c(length(rnaTypes),5))
  br = seq(0, maxlen, by=200)
  for(rt in rnaTypes) {
    s   = cs[[rt]]$s
    stopifnot(all(selectedCategories %in% levels(s$category)))
    for(lev in selectedCategories) {
      len = s$length[s$category == lev]
      len[len>maxlen]=maxlen
      hist(len, breaks=br, col=fillColors[lev], main=paste(rt, lev))
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
    pdf("tableSegments-lvsx.pdf", width=14, height=length(rnaTypes)*3)
    pch="."
  } else {
    pch=18
  }
  par(mfrow=c(length(rnaTypes), 2))
  maxlen=5000
  br = seq(0, maxlen, by=200)
  selectedCategories = c("verified", "unIso")
  for(rt in rnaTypes) {
    s   = cs[[rt]]$s
    stopifnot(all(selectedCategories %in% levels(s$category)))
    ylim = quantile(s$level[s$category %in% selectedCategories], probs=c(0.01, 0.99), na.rm=TRUE)
    for(lev in selectedCategories) {
      len = s$length[s$category == lev]
      len[len>maxlen] = maxlen
      exl = s$level[s$category == lev]
      
      plot(len, exl, pch=pch, ylim=ylim, main=paste(longNames[rt], ": ", lev, sep=""),
           ylab="expression level", xlab="length")
      lf = loess(exl ~ len)
      slen = sort(len)
      lines(slen, predict(lf, newdata=slen), col="blue")
    }
    rm(list=c("s", "ylim", "len", "lev", "exl", "lf", "slen"))
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
    
calchit = function(sp, blrt, s) {  
  hit = matrix(NA, nrow=length(sp), ncol=length(blrt))
  rownames(hit) = names(sp)
  colnames(hit) = names(blrt)
  for(b in 1:length(blrt)) {
    br  = blrt[[b]]
    ## 1 = Query Sequence ID, 3 = Percent identity, 4 = Alignment length 
    fas = br[[3]] * br[[4]] / (s$length[br[[1]]]+1)
    ## fas = br[[4]] / (s$length[br[[1]]]+1) * 100
    ## fas = br[[3]] 
    stopifnot(all( fas>=0 & fas<=115 & !is.na(fas)))
    
    ## split by name of query sequence and just keep the hit with the
    ## highest value of fas
    spbyq   = split(1:nrow(br), br[[1]]) 
    theBest = sapply(spbyq, function(i) i[which.max(fas[i])])
    i.seg = br[[1]][theBest]
      
    ## mean of the ratios:
    alignableFrac = numeric(nrow(s))
    alignableFrac[i.seg] = fas[theBest] 
    
    hit[,b] = sapply(sp, function(segments) {
      mean(alignableFrac[segments]) 
    })
  }
  data.frame(number=listLen(sp), fraction=rowMeans(hit))
}

##
## cons
##

if("cons" %in% what){
  cat("Fraction of alignable sequence (percent):\n",
      "=========================================\n",
      "Here, 'number' is the number of segments.\n", sep="")

  theSplit = vector(mode="list", length=length(rnaTypes))
  names(theSplit) = rnaTypes
    
  selectedCategories = c("verified", "uncharacterized" , "ncRNA",  "unAnti", "unIso")
  
  for(rt in rnaTypes) {
    s   = cs[[rt]]$s
    stopifnot(all(selectedCategories %in% levels(s$category)))
    sp = split(seq(along=s$category), s$category)
    sp = sp[selectedCategories]
    sp$"unannotated, unexpressed"  = which(s$overlappingFeature=="" &
        s$oppositeFeature=="" & s$level <= quantile(s$level, 0.2, na.rm=TRUE))
    sp$"whole genome" = seq(along=s$category)
    
    hit = calchit(sp, blastres[[rt]], s)

    theSplit[[rt]] = sp
    cat("\n", rt, "\n-----\n", sep="")
    print(round(hit,1))
  }

  nrlevs = 4
  minlen = 200
  cat("\n\nGrouping of conservation scores by expression, using ", nrlevs, "\n", 
      "quantile groups of equal size, respectively.\n",
      "Scores are calculated for segments with length >= ", minlen, ".\n",
      "==========================================================================\n", 
      sep="")
  
  for(plotWhat in 1:2) {
    if(interact) {
      x11(width=6*length(rnaTypes), height=7)
    } else {
      pdf(paste("tableSegments-consex-", plotWhat, ".pdf", sep=""),
          width=4*length(rnaTypes), height=7)
    }
    par(mfcol=c(1, length(rnaTypes)))
    
    for(rt in rnaTypes) {
      cat("\n", rt, "\n-----\n", sep="")

      sp = switch(plotWhat,
        theSplit[[rt]][c("verified", "unIso")],
        theSplit[[rt]]
        )
        
      names(longNames) = longNames = names(sp)
      longNames[names(longNames)=="verified"] = "verified genes"
      longNames[names(longNames)=="unIso"]    = "unannotated, isolated"
                
      fraction = matrix(NA, nrow=nrlevs, ncol=length(sp))
      colnames(fraction) = names(sp)
  
      pchs = switch(plotWhat,
        14+seq(along=sp),
        rep(16, length(sp))
        )
      
      s   = cs[[rt]]$s
      for(j in seq(along=sp)) {
        wh   = sp[[j]][s$length[sp[[j]]]>=minlen]
        v    = s$level[wh]
        br   = quantile(v, (0:nrlevs)/nrlevs, na.rm=TRUE)
        spwh = split(wh, cut(v, breaks=br))
        hit = calchit(spwh, blastres[[rt]], s)
        fraction[, names(sp)[j]] = hit$fraction
        cat("\nby expression for:", names(sp)[j], "\n")
        print(round(hit,1))
      }
      matplot(fraction, xaxt="n", type="b", main=rt,
              lty=1, lwd=2, pch=pchs, col=lineColors[colnames(fraction)],
              ylab="average identity (percent) ", xlab="transcript level",
              ylim=c(0, 110))
      cutNames = paste("<=", round((1:nrlevs)/nrlevs*100), "%", sep="")
      axis(side=1, at = 1:nrlevs, labels = cutNames)
      if(rt=="polyA2")
        legend(x=1, y=112, legend=longNames[colnames(fraction)],
               lty=1, lwd=2, pch=pchs, col=lineColors[colnames(fraction)])
    }
    if(!interact)
      dev.off()
  }
}

if(!interact)
  sink()
