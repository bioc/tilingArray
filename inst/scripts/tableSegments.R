library("tilingArray")
library("geneplotter"); source("~/madman/Rpacks/geneplotter/R/histStack.R")

graphics.off()
options(error=recover, warn=2)
interact = (!TRUE)
what     = c("fig2", "wpt", "cons", "lvsx", "wst")[1:4] 

consScoreFun = function(alignmentLength, percentIdentity, queryLength)
  (alignmentLength*percentIdentity/queryLength)

rnaTypes  = c("seg-polyA-050525", "seg-tot-050525", "seg-tot2-050525")[1:2]
outfile = "tableSegments"

source("scripts/readSegments.R") 
source("scripts/categorizeSegments.R") 
source("scripts/writeSegmentTable.R")

if(!interact){
  sink(paste(outfile, ".txt", sep=""))
  cat("Made on", date(), "\n\n")
}

source("scripts/calcThreshold.R") 

##
## CATEGORIZE
##
if(!exists("cs")) {
  cs = vector(mode="list", length=length(rnaTypes))
  names(cs)=rnaTypes

  cat("Categorization of segments:\n",
      "===========================\n", sep="")
  for(rt in rnaTypes) {
    cat(rt, ":\n", sep="")
    s = categorizeSegments(get(rt))
    cs[[rt]] =s
  }
} else {
  cat("\n**************************************************\n",
        "*      NOT REDOING categorizeSegments            *\n",
        "**************************************************\n", sep="")
}
  
fillColors = c(c(brewer.pal(10, "Paired")[c(1, 2, 6, 8, 5:8, 2, 10)]), "#d0d0d0")
names(fillColors) = c("overlap <50%", "overlap >=50%",
       "novel antisense", "novel isolated",
       "novel isolated - unassigned",  "novel isolated - filtered",
       "novel antisense - unassigned", "novel antisense - filtered",
         "annotated ORF", "ncRNA(all)", "untranscribed")

lineColors = c(brewer.pal(8, "Paired")[c(1,2,6,8)], "grey")
names(lineColors) =c("annotated ORF", "ncRNA(all)", 
       "novel antisense - filtered", "novel isolated - filtered",
       "unexpressed isolated")

##
## PIE: Four classes
##
if("fig2" %in% what){

  if(interact) {
    x11(width=4*length(rnaTypes), height=3*4)
  } else {
    pdf(paste(outfile, "fig2.pdf", sep="-"), width=4*length(rnaTypes), height=3*4.2)
  }

  par(mfrow=c(3, length(rnaTypes)))
  counts = NULL

  cat("\nSegment overlap with known features (genes):\n",
        "============================================\n\n", sep="")

  pieCat = vector(mode="list", length=length(rnaTypes))
  names(pieCat) = rnaTypes
  
  for(irt in seq(along=rnaTypes)) {
    rt = rnaTypes[irt]
    s  = cs[[rt]] 
   
    px = table(s[, "pieCat"])
    labels = LETTERS[ match(names(px), levels(s[, "pieCat"])) ]
    stopifnot(all(names(px) %in% names(fillColors)))
    counts = cbind(counts, px)
    pie(px, radius=0.75, main=paste(c("a)", "b)")[irt], longNames[rt]), 
        col = fillColors[names(px)], labels = paste(labels, " (", px, ")", sep=""))

    category = s[, "category"]
    levels(category) = sub("ncRNA", "ncRNA(all)", levels(category))
    category[ s[, "simpleCatg"]=="ncRNA(all)" ] = "ncRNA(all)"

    cat(rt, ":\n", sep="")
    tab = table(category, s[, "overlap"])
    tab = tab[rowSums(tab)!=0, ]
    print(tab)
    cat("\n\n")
  } ## for rt

  colnames(counts)=rnaTypes
  cat("\nSegment counts (pie charts):\n",
        "============================\n\n", sep="")
  print(counts)
  cat("\n")
  
  ##
  ## Compare total to poly-A, the goal is: which transcripts do we find 
  ## specifically in total RNA?
  ##
  s1     = cs[[rnaTypes[1]]]
  s2     = cs[[rnaTypes[2]]]
  start1 = s1[, "start"]
  end1   = s1[, "end"]
  start2 = s2[, "start"]
  end2   = s2[, "end"]

  unTrCatgs = c("excluded", "untranscribed")
  ## unTrCatgs = c(unTrCatgs, "novel isolated - unassigned",
  ##  "novel antisense - unassigned")
    
  isTr1  = !(s1[, "category"] %in% unTrCatgs)
  
  jstart = jend = 1
  ov     = numeric(nrow(s2))
  for(k in 1:nrow(s2)) {
    ## make sure that jstart points to a segment in s1 whose
    ## start <= ks <= end, where ks=start of current segment in s2.
    ks = start2[k]
    ke = end2[k]
    while(!((start1[jstart]<=ks) && (end1[jstart]>=ks)))
      jstart = jstart+1
    ## Similarly, make sure that jend points to a segment in s1 whose
    ## start <= ke <= end, where ke=end of current segment in s2.
    while(!((start1[jend]<=ke) && (end1[jend]>=ke)))
      jend = jend+1
    stopifnot(jstart<=nrow(s1), jend<=nrow(s1))
    ## cat(k, jstart, jend, "\n")
    
    if(jstart==jend) {
      lne = lns = 0
      lni = ke-ks+1
    } else {
      lns = (end1[jstart] - ks + 1) * as.numeric(isTr1[jstart])
      lne = (ke - start1[jend] + 1) * as.numeric(isTr1[jend])
      if(jend-jstart>1) {
        j = (jstart+1):(jend-1)
        lni = sum( (end1[j]-start1[j]+1) * as.numeric(isTr1[j]) )
      } else {
        lni = 0
      }
    }
    stopifnot(lns>=0, lne>=0, lni>=0)
    ov[k] = (lns+lni+lne) / (ke-ks+1)
  }
  
  ##if(!interact)
  ##  pdf(paste(outfile, "overlap.pdf", sep="-"), width=7, height=4.8)
  ##hist(ov, 100, col="orange", xlab="overlap", main="")
  ##if(!interact)
  ##  dev.off()
  
  cat("\n\nOverlap of segments from total RNA with those from poly-A\n",
          "(other than:", paste(unTrCatgs, collapse=", "), ")\n",
          "=========================================================\n", sep="")
  
  tab = table(s2[, "category"], ov > .5)
  tab = tab[ !(rownames(tab)%in%c("excluded", "untranscribed")), 2:1]
  print(tab)
  cat("\n\n")

  ##
  ##
  selectedCategories = c(
    "(1) only 'overlap <50%' (i.e. in 'overlappingFeature' but not 'mostOfFeatureinSegment')", 
    "(2) only 'overlap >=50%', but not 'complete' (i.e. in 'mostOfFeatureinSegment' but not 'featureInSegment')",
    "(3) 'complete' (i.e. in 'featureInSegment')",
    "(4): (1) or (2)",
    "(5): (2) or (3)")

  selGenes = ((gff[,"feature"]=="gene") & (gff[, "orf_classification"] %in% c("Verified", "Uncharacterized")))
  featNames = list("annotated ORFs" = unique(gff[ selGenes, "Name"]),
                   "ncRNA(all)" = unique(gff[ gff[, "feature"] %in% allncRNA, "Name"]))

  tab = matrix(NA, nrow=length(selectedCategories)*length(featNames), ncol=length(rnaTypes)+1)
  rownames(tab) = paste("(", rep(1:5, 2), ") ", rep(names(featNames), each=5), sep="")
  colnames(tab) = c(rnaTypes, "in genome")

  for(rt in rnaTypes) {
    ovf = unique(unlist(strsplit(cs[[rt]][, "overlappingFeature"],     split=", ")))
    mof = unique(unlist(strsplit(cs[[rt]][, "mostOfFeatureInSegment"], split=", ")))
    fis = unique(unlist(strsplit(cs[[rt]][, "featureInSegment"],       split=", ")))
    
    for(isc in seq(along=selectedCategories)) {
      fIDs = switch(isc,
        setdiff(ovf, mof),
        setdiff(mof, fis),
        fis,
        setdiff(ovf, fis),
        mof)
      for(k in seq(along=featNames))
        tab[ (k-1)*length(selectedCategories)+isc, rt ] = length(intersect(fIDs, featNames[[k]]))
    }    
  }
  tab[ , "in genome" ] = rep(listLen(featNames), each=length(selectedCategories))

  cat("\nHow many unique known features (SGD Names) do we find with\n",
      paste(selectedCategories, collapse="\n"), "\n", 
      "=====================================================================\n", sep="")
  print(tab)  
  cat("\n\n")

  ##
  ## LENGTH & LEVEL DISTRIBUTIONS
  ##
  maxlen=5000
  br = seq(0, maxlen, by=200)
  breaksFun = function(z) paste(signif(z, 3))
  for(irt in seq(along=rnaTypes)) {
    s   = cs[[rnaTypes[irt]]]
    plotCat = s[, "pieCat"]
    stopifnot(all(levels(plotCat) %in% names(fillColors)))
    len = split(s[, "length"], plotCat)
    len = lapply(len, function(z) {z[z>maxlen]=maxlen; z})
    len = len[order(listLen(len))]
    cols = fillColors[names(len)]
    histStack(len, breaks=br, col=cols, 
              main=paste(c("c)", "d)")[irt]), breaksFun=breaksFun, xlab="length")
  }

  for(irt in seq(along=rnaTypes)) {
    s  = cs[[rnaTypes[irt]]]
    plotCat = s[, "pieCat"]
    levels(plotCat) = c(levels(plotCat), "untranscribed")
    plotCat[ s[, "simpleCatg"]=="untranscribed" ] = "untranscribed"
    stopifnot(all(levels(plotCat) %in% names(fillColors)))

    thr = get("threshold", get(rnaTypes[irt]))
    lv = split(s[, "level"], plotCat)
    lv = lv[order(listLen(lv))]
    by = 0.2
    rg = quantile(unlist(lv), probs=c(0.001, 0.999))
    br = c(rev(seq(thr, rg[1]-by, by=-by)), seq(thr+by, rg[2]+by, by=by))
    lv = lapply(lv, function(z)
      ifelse(z<=rg[2], ifelse(z>=rg[1], z, rg[1]), rg[2]))
    histStack(lv, breaks=br, col=fillColors[names(lv)],
              main=paste(c("e)", "f)")[irt]), breaksFun=breaksFun, xlab="level")
    ## abline(v=which.min(abs(br-thr))-1, col="grey", lwd=3)
  }
  
  if(!interact)
    dev.off()
}

##
## What fraction of basepairs in the genome are transcribed
## and how many genes do we find expressed
##
if("wpt" %in% what){
  data(transcribedFeatures)
  feats = transcribedFeatures[transcribedFeatures!="gene"]
  nrChr = 16

  chrlen = sapply(1:nrChr, function(chr)
    max(gff[gff[, "chr"]==chr , "end"]))

  isAnno = lapply(1:nrChr, function(chr) {
    res  = logical(chrlen[chr])
    selg = which((gff[, "chr"]==chr) & (gff[, "feature"] %in% transcribedFeatures))
    for(j in selg)
      res[gff$start[j]:gff$end[j]] = TRUE
    res
  })
  isAnno = unlist(isAnno)

  isTrans = segLev = vector(mode="list", length=3)
  names(isTrans) = names(segLev) = c(rnaTypes, "both")
  
  for(rt in rnaTypes) {
    s = cs[[rt]]
    threshold = get("threshold", get(rt))
    res = lapply(1:nrChr, function(chr) {
      res  = rep(-Inf, chrlen[chr])
      selt = which(s[, "chr"]==chr & !is.na(s[, "level"]) & s[,"frac.dup"]<maxDuplicated)
      for(i in selt) {
        rg = s$start[i]:s$end[i] 
        res[rg] = pmax(res[rg], s[i, "level"])
      }
      res
    })
    segLev[[rt]]  = unlist(res)
    isTrans[[rt]] = (segLev[[rt]] >= get("threshold", get(rt)))
  }
  isTrans[["both"]] = (isTrans[[1]] | isTrans[[2]])
  segLev[["both"]]  = pmax(segLev[[1]], segLev[[2]])
  
  cat("Fraction of transcribed basepairs\n",
      "=================================\n\n", sep="")
  percent ="%"
  cat(sprintf("%31s: %7d of %7d bp (%3.1f%s)\n\n", "Annotated", sum(isAnno), length(isAnno),
      signif(mean(isAnno)*100, 3), percent))
  for(i in seq(along=isTrans)) {
    n1    = sum(isTrans[[i]])
    n2    = sum(isTrans[[i]] & !isAnno)
    denom = sum(is.finite(segLev[[i]]))
    cat(sprintf("Transcribed in %16s: %7d of %7d bp (%3.1f%s)\n", names(isTrans)[i], 
      n1, denom, signif(n1/denom*100, 3), percent))
    cat(sprintf("          ... and not annotated: %7d of %7d bp (%3.1f%s)\n\n", 
      n2, denom, signif(n2/denom*100, 3), percent))
    
  }
  cat("\n")

  if(!interact)
    pdf(paste(outfile, "wpt.pdf", sep="-"), height=length(rnaTypes)*4, width=6)
  par(mfrow=c(2,1))
  cols = brewer.pal(4, "Paired")
  for(i in seq(along=rnaTypes)) {
    rt = rnaTypes[i]
    hist(segLev[[rt]], breaks=100, col=cols[i*2-1], main=rt, xlab="level")
    abline(v=get("threshold", get(rt)), col=cols[i*2], lwd=3)
  }
  if(!interact)
    dev.off()  
}

##
## WRITE THE SEGMENT TABLE
##
if("wst" %in% what){
  for(rt in rnaTypes) {
    s = cs[[rt]]
    s$segID = paste(1:nrow(s))
    drop =  (s[,"category"]=="excluded") | (s[,"category"]=="untranscribed"&(!s[,"isUnIso"]))
    writeSegmentTable(s[!drop, ],
      fn = file.path(indir[rt], "viz", "index"), HTML=TRUE, 
      sortBy = "category-level",
      title = paste(rt, " (", longNames[rt], ")", sep=""), interact=interact)
  }
  cat("\n")
}

##
## LENGTH VERSUS EXPRESSION LEVEL
##
if("lvsx" %in% what){
  if(!interact) {
    pdf(paste(outfile, "lvsx.pdf", sep="-"), width=14, height=length(rnaTypes)*3)
    pch="."
  } else {
    pch=18
  }
  par(mfrow=c(length(rnaTypes), 2))
  maxlen=5000
  br = seq(0, maxlen, by=200)
  selectedCategories = c("annotated ORF", "novel isolated - filtered")
  for(rt in rnaTypes) {
    s = cs[[rt]]
    stopifnot(all(selectedCategories %in% levels(s[,"simpleCatg"])))
    ylim = quantile(s[,"level"][s[,"simpleCatg"] %in% selectedCategories],
      probs=c(0.01, 0.99), na.rm=TRUE)
    for(lev in selectedCategories) {
      len = s[s[,"simpleCatg"] == lev, "length"]
      exl = s[s[,"simpleCatg"] == lev, "level"]
      len[len>maxlen] = maxlen
      plot(len, exl, pch=pch, ylim=ylim,
           main=paste(longNames[rt], ": ", lev, sep=""),
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
if("cons" %in% what){

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
    fas = consScoreFun(br[[4]], br[[3]], s$length[br[[1]]])
    ##fas = (br[[4]] / s$length[br[[1]]] > 0.5) * 100 
    ##fas = rep(100, nrow(br))
    stopifnot(all( fas>=0 & fas<=115 & !is.na(fas)))
    
    ## split by name of query sequence and just keep the hit with the
    ## highest value of fas
    spbyq   = split(1:nrow(br), br[[1]]) 
    theBest = sapply(spbyq, function(i) i[which.max(fas[i])])
    i.seg = br[[1]][theBest]
    
    ## mean of the ratios:
    alignableFrac = numeric(nrow(s))
    alignableFrac[i.seg] = fas[theBest]

    ## number of hits somewhere else - this seems to be non-sense
    ## tab = table(br[[1]])
    ## alignableFrac = numeric(nrow(s))
    ## alignableFrac[as.numeric(names(tab))]=tab
    
    hit[,b] = sapply(sp, function(segments) {
      mean(alignableFrac[segments]) 
    })
  }
  cbind(number=listLen(sp), hit)
}

  cat("Fraction of alignable sequence (percent):\n",
      "=========================================\n",
      "Here, 'number' is the number of segments.\n", sep="")

  theSplit = vector(mode="list", length=length(rnaTypes))
  names(theSplit) = rnaTypes
    
  selectedCategories = c("annotated ORF", "ncRNA(all)",
    "novel antisense - filtered", "novel isolated - filtered",
    "unexpressed isolated")

  for(rt in rnaTypes) {
    s   = cs[[rt]]

    catg = s[,"simpleCatg"]
    levels(catg) = c(levels(catg), "unexpressed isolated")
    catg[s[,"isUnIso"]] = "unexpressed isolated"

    stopifnot(all(selectedCategories %in% levels(catg)))
    sp = split(seq(along=catg), catg)
    sp = sp[selectedCategories]
    
    hit = calchit(sp, blastres[[rt]], s)
    theSplit[[rt]] = sp
    cat("\n", rt, "\n-----\n", sep="")
    print(round(hit,1))
  }

  nrlevs = 3
  cat("\n\nGrouping of conservation scores by expression,\nusing ", nrlevs, 
      " quantile groups of equal size.\n",
      "=====================================================================\n", 
      sep="")

  species = names(blastres[[1]])
  if(interact) {
    x11(width=4.5*length(rnaTypes), height=4.5*length(species))
  } else {
    pdf(paste(outfile, "consex.pdf", sep="-"),
        width=4.5*length(rnaTypes), height=4.5*length(species))
  }
  par(mfcol=c(length(species), length(rnaTypes)))
  stopifnot(all(selectedCategories %in% names(lineColors)))

  for(i in seq(along=rnaTypes)) {
    rt=rnaTypes[i]
    cat("\n", rt, "\n------------\n", sep="")
    
    sp = theSplit[[rt]]
    s  = cs[[rt]]
    
    fraction = array(NA, dim=c(nrlevs, length(species), length(sp)))
    dimnames(fraction)=list(NULL, species=species, category=names(sp))
  
    pchs = c(15, 5, 16, 17, 19) # 20
    stopifnot(length(pchs)==length(sp))
    
    for(j in seq(along=sp)) {
      wh   = sp[[j]]
      v    = s$level[wh]
      br   = quantile(v, (0:nrlevs)/nrlevs, na.rm=TRUE)
      spwh = split(wh, cut(v, breaks=br))
      hit  = calchit(spwh, blastres[[rt]], s)
      fraction[,, names(sp)[j]] = hit[, species]
      cat("\nby expression for:", names(sp)[j], "\n")
      print(round(hit,1))
    }
    #j = which(colnames(fraction)=="isolated and unexpressed")
    #stopifnot(length(j)==1)
    #baseline = mean(fraction[, j])
    #fraction = fraction[, -j]

    for(j in seq(along=species)){
      linNams = dimnames(fraction)[[3]]
      linCols = lineColors[ linNams ]
      matplot(fraction[,j,], xaxt="n", type="b", main=paste(species[j], "--", longNames[rt]),
            lty=1, lwd=2, pch=pchs, col=linCols,
            ylab="average identity (percent) ", xlab="transcript level",
            ylim=c(0, 110))
    
      ## abline(h=baseline, col=lineColors[["isolated and unexpressed"]], lwd=2, type="l")
    
      cutNames = paste("<=", round((1:nrlevs)/nrlevs*100), "%", sep="")
      axis(side=1, at = 1:nrlevs, labels = cutNames)
      if(i==1)
        legend(x=1, y=112, legend=linNams, lty=1, lwd=2, pch=pchs, col=linCols)
    } ## for j
  } ## for i
  if(!interact)
    dev.off()
 
}

if(!interact)
  sink()
