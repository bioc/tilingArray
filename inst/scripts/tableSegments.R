library("tilingArray")
source("scripts/readSegments.R") 
source("scripts/calcThreshold.R") 
source("scripts/categorizeSegments.R") 
source("scripts/writeSegmentTable.R")

options(error=recover, warn=2)
## debug(categorizeSegmentsPie)

interact =  FALSE
what=c("pie", "wst", "length", "lvsx", "cons", "conswex")

if(!interact & exists("cs"))
  rm(cs)

n = length(rnaTypes)

if(!interact)
  sink("tableSegments.txt")

if(!exists("cs")) {
  graphics.off()
  if(interact)
    x11(width=10, height=n*3)
  else
    pdf(width=11, height=n*4)
  par(mfrow=c(n,1))
  
  cs = vector(mode="list", length=length(rnaTypes))
  names(cs)=rnaTypes
  
  cat("Categorizing segments:\n",
      "======================\n", sep="")
  for(rt in rnaTypes) {
    s = get("segScore", get(rt))
    cs[[rt]] = categorizeSegmentsPie(s)
  }

  if(!interact)
    dev.off()
}


colors = c(brewer.pal(9, "Pastel1")[c(2:6, 1)])
names(colors) =c("verified", "ncRNA", "uncharacterized", "dubious",
                 "unAnti", "unIso")

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
    cat("\n", rt, "\n-----\n", sep="")
    selectedCategories = c("verified", "ncRNA", "uncharacterized", "dubious", "unIso", "unAnti")
    stopifnot(all(selectedCategories %in% names(colors)))
    ct  = cs[[rt]]$count
    stopifnot(all(selectedCategories %in% rownames(ct)))

    print(cbind(ct, "percent" = round(100*ct[,1]/ct[,2], 1)))

    px = ct[selectedCategories, "observed"]
    pie(px, radius=0.75, main=longNames[rt],
        col=colors[selectedCategories],
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
  stopifnot(all(selectedCategories %in% names(colors)))
  if(!interact)
    pdf("tableSegments-lengths.pdf", width=14, height=n*3)
  par(mfrow=c(n,5))
  br = seq(0, maxlen, by=200)
  for(rt in rnaTypes) {
    s   = cs[[rt]]$s
    stopifnot(all(selectedCategories %in% levels(s$category)))
    for(lev in selectedCategories) {
      len = s$length[s$category == lev]
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
  cat("\n\nFraction of alignable sequence (percent):\n",
      "=========================================\n",
      "Here, 'number' is the number of segments.\n", sep="")
  
  selectedCategories = c("verified", "uncharacterized" , "ncRNA",  "unAnti", "unIso")
  for(rt in rnaTypes) {
    s   = cs[[rt]]$s
    stopifnot(all(selectedCategories %in% levels(s$category)))
    sp = split(seq(along=s$category), s$category)
    sp = sp[selectedCategories]
    sp$"unexpressed & unannotated"  = which(s$overlappingFeature=="" &
        s$oppositeFeature=="" & s$level <= quantile(s$level, 0.2, na.rm=TRUE))
    sp$"whole genome" = seq(along=s$category)
    
    hit = calchit(sp, blastres[[rt]], s)
    
    cat("\n", rt, "\n-----\n", sep="")
    print(round(hit,1))
  }
}

##
## conswex
##
if("conswex" %in% what){
  nrlevs = 3
  minlen = 150
  evList = c("expression", "length")[1]
  cat("\n\nGrouping of conservation scores by ", paste(evList, collapse=", "), "\n",
      "using ", nrlevs,
      " quantile groups of equal size, respectively.\n",
      "For 'expression', scores are calculated for segments with length >= ", minlen, ".\n",
      "==========================================================================\n", 
      sep="")
  if(!interact)
    pdf("tableSegments-conswex.pdf", width=4*n, height=4*length(evList))
  par(mfcol=c(length(evList), n))
  cols = c("#303030", "#0000e0")
  pchs = c(15,16)
  
  selectedCategories = c("verified", "unIso")
  selectedCategoriesLong = c("verified genes", "unannotated, isolated")
  fraction = matrix(NA, nrow=nrlevs, ncol=length(selectedCategories))
  colnames(fraction) = selectedCategories
  
  for(rt in rnaTypes) {
    s   = cs[[rt]]$s
    stopifnot(all(selectedCategories %in% levels(s$category)))
    cat("\n", rt, "\n-----\n", sep="")

    for(ev in evList) {
      for(segClass in selectedCategories) {
        wh = which(s$category == segClass)
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

        hit = calchit(sp, blastres[[rt]], s)
        fraction[, segClass] = hit$fraction
        cat("\nby '", ev, "', for '", segClass, "' segments\n", sep="")
        print(round(hit,1))
      }
      matplot(fraction, xaxt="n", type="b", main=paste(ev, " (", rt, ")", sep=""),
              lty=1, lwd=2, pch=pchs, col=cols,
              ylab="Fraction of alignable sequence", xlab=ev,
              ylim=c(0, 90))
      axis(side=1, at = 1:nrlevs, labels = names(sp))
      if(rt=="polyA2" & ev=="expression")
        legend(x=1, y=90, legend=selectedCategoriesLong, lty=1, lwd=2, pch=pchs, col=cols)
    }
  }
  if(!interact)
    dev.off()
}

if(!interact)
  sink()
