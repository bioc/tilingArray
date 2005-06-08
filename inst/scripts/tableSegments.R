library("tilingArray")

graphics.off()
options(error=recover, warn=2)
interact = (TRUE)
what     = c("pie", "wpt", "wst", "length", "cons", "lvsx")[3] 

consScoreFun = function(alignmentLength, percentIdentity, queryLength)
  (alignmentLength*percentIdentity/queryLength)

source("scripts/readSegments.R") 
source("scripts/categorizeSegments.R") 
source("scripts/writeSegmentTable.R")

if(!interact){
  outfile = "tableSegments"
  ## outfile = "tableSegments-yesno"
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
  
  for(rt in rnaTypes) {
    cat("\n--------", rt, "---------\n")
    s = categorizeSegments(get(rt))
    cs[[rt]] =s
  }
  
  fillColors = c(brewer.pal(10, "Paired")[c(1, 2, 6, 8, 5:8, 2, 10)])
  names(fillColors) = c("overlap < 50%", "overlap >=50%", "novel antisense", "novel isolated",
         "novel antisense - excluded", "novel antisense - filtered",
         "novel isolated - excluded", "novel isolated - filtered",
         "annotated ORF", "ncRNA(all)")
  
  lineColors = c(brewer.pal(8, "Paired")[c(1,2,6,8)], "grey")
  names(lineColors) =c("annotated ORF", "ncRNA(all)", 
   "novel antisense - filtered", "novel isolated - filtered", "unexpressed isolated")
} else {
  cat("\n**************************************************\n",
        "*      NOT REDOING categorizeSegments            *\n",
        "**************************************************\n", sep="")
}
##
## PIE: Four classes
##
if("pie" %in% what){

  if(interact) {
    x11(width=7*length(rnaTypes), height=4.8)
  } else {
    pdf(paste(outfile, "pie.pdf", sep="-"), width=7*length(rnaTypes), height=4.8)
  }

  par(mfrow=c(1, length(rnaTypes)))
  counts = NULL
  cat("\n\n")
  for(rt in rnaTypes) {
    s  = cs[[rt]] 
    category = s[, "category"]
    overlap  = s[, "overlap"]

    px = c("overlap >=50%"         = sum(overlap  %in% c(">=50%, <100%", "100%")),
      "overlap < 50%"              = sum(overlap  %in% c("<50%")),
      "novel isolated - filtered"  = sum(category %in% c("novel isolated - filtered")),
      "novel isolated - excluded"  = sum(category %in% c("novel isolated - excluded")),
      "novel antisense - filtered" = sum(category %in% c("novel antisense - filtered")),
      "novel antisense - excluded" = sum(category %in% c("novel antisense - excluded")))
##   c("overlap >=50%"   = sum(overlap  %in% c(">=50%, <100%", "100%")),
##     "overlap < 50%"   = sum(overlap  %in% c("<50%")),
##     "novel isolated"  = sum(category %in% c("novel isolated - filtered", "novel isolated - excluded")),
##     "novel antisense" = sum(category %in% c("novel antisense - filtered", "novel antisense - excluded"))),
    
    stopifnot(all(names(px) %in% names(fillColors)))
    counts = cbind(counts, px)
    pie(px, radius=0.75, main=longNames[rt], 
        col = fillColors[names(px)],
        labels = paste(names(px), " (", px, ")", sep=""))

    cat("=====", rt, "=====\n")
    levels(category) = sub("ncRNA", "ncRNA(all)", levels(category))
    category[ s[, "simpleCatg"]=="ncRNA(all)" ] = "ncRNA(all)"
    tab = table(category, overlap)
    tab = tab[rowSums(tab)!=0, ]

    ## add a column "in genome"
    featg = as.character(gff[, "feature"])
    for(grx in c("Dubious", "Verified", "Uncharacterized"))
      featg[ (featg=="gene") & (gff[, "orf_classification"]==grx) ] = paste(tolower(grx), "gene")
    featg[ featg %in% allncRNA ] = "ncRNA(all)"
    
    tab = cbind(tab, "in genome" = sapply(rownames(tab), function(ft)
                       length(unique(gff[ featg==ft, "Name"]))))
      
    print(tab)
    cat("\n\n")
  }
  if(!interact)
    dev.off()

  colnames(counts)=rnaTypes
  rownames(counts)=names(px)
  cat("\nSegment counts:\n")
  print(counts)
  cat("\n\n")
} ## if


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

  segLev = isTrans = isUni = vector(mode="list", length=length(rnaTypes))
  names(segLev) = names(isTrans) = names(isUni) = rnaTypes
  
  for(rt in rnaTypes) {

    s = cs[[rt]]
    threshold = get("threshold", get(rt))
      
    res = lapply(1:nrChr, function(chr) {
      res  = rep(-Inf, chrlen[chr])
      selt = which(s[, "chr"]==chr & !is.na(s[, "level"])) 
      for(i in selt) {
        rg = s$start[i]:s$end[i] 
        res[rg] = pmax(res[rg], s[i, "level"])
      }
      res
    })
    segLev[[rt]] = unlist(res)
    isTrans[[rt]] = (segLev[[rt]] >= get("threshold", get(rt)))

    res = lapply(1:nrChr, function(chr) {
      res  = logical(chrlen[chr])
      selt = which(s[, "chr"]==chr & s$frac.dup>=maxDuplicated)
      for(i in selt)
        res[s$start[i]:s$end[i]] = TRUE
      res
    })
    isUni[[rt]] = !(unlist(res))
  }
  
  cat("Fraction of transcribed basepairs\n",
      "=================================\n\n", sep="")
  cat(sprintf("%31s: %3.1f percent\n", "Annotated", signif(mean(isAnno)*100, 3)))
  for(rt in rnaTypes) {
    cat(sprintf("Transcribed in %16s: %3.1f percent\n", rt, 
     signif(sum(isTrans[[rt]])/sum(isUni[[rt]])*100, 3)))
  }
  cat(sprintf("%31s: %3.1f percent\n", "Union of both",
     signif(sum(isTrans[[1]]|isTrans[[2]])/sum(isUni[[1]]|isUni[[2]])*100, 3)))
  cat("\n\n")

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
## LENGTH & LEVEL DISTRIBUTIONS
##
if("length" %in% what){
  selectedCategories = c(
     "annotated ORF", "ncRNA(all)", 
     "novel isolated - filtered", 
     "novel antisense - filtered")

  stopifnot(all(selectedCategories %in% names(fillColors)))
  
  if(!interact)
    pdf(paste(outfile, "lengths.pdf", sep="-"), width=14, height=length(rnaTypes)*3)
  par(mfrow=c(length(rnaTypes), length(selectedCategories)))
  maxlen=5000
  br = seq(0, maxlen, by=200)
  for(rt in rnaTypes) {
    s = cs[[rt]]
    stopifnot(all(selectedCategories %in% levels(s[,"simpleCatg"])))
    for(lev in selectedCategories) {
      len = s[ s[, "simpleCatg"]==lev, "length"]
      len[len>maxlen]=maxlen
      hist(len, breaks=br, col=fillColors[lev], main=paste("Length:", rt, lev), xlab="length")
    }
  }
  if(!interact) {
    dev.off()
    pdf(paste(outfile, "levels.pdf", sep="-"), width=14, height=length(rnaTypes)*3)
  }
  par(mfrow=c(length(rnaTypes), length(selectedCategories)))
  br = seq(get("threshold", get(rt))-0.01, max(s[, "level"], na.rm=TRUE)+0.01, length=40)
  for(rt in rnaTypes) {
    s = cs[[rt]]
    for(lev in selectedCategories) {
      lvs = s[ s[, "simpleCatg"]==lev, "level"]
      hist(lvs, breaks=br, col=fillColors[lev], main=paste("Level:", rt, lev), xlab="level")
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

##
## cons
##

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
