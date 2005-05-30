library("tilingArray")

graphics.off()
options(error=recover, warn=2)
interact = (!TRUE)
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

for(rt in rnaTypes) {
  cs[[rt]] = categorizeSegments(get(rt))
}

#fillColors = , 1, 9)])
#names(fillColors) =c("annotated ORF", "ncRNA", "other annotation",
#   "novel antisense", "novel isolated", "excluded")

fillColors = c(brewer.pal(8, "Paired")[c(1,2,6,8)], brewer.pal(8, "Paired")[5:8],
               brewer.pal(9, "Pastel1")[c(2:3)])
names(fillColors) = c("overlap < 50%", "overlap >=50%", "novel antisense", "novel isolated",
       "novel antisense - excluded", "novel antisense - filtered",
       "novel isolated - excluded", "novel isolated - filtered",
       "annotated ORF", "ncRNA")
      
lineColors = c(brewer.pal(8, "Paired")[c(1,2,6,8)], "grey")
names(lineColors) =c("annotated ORF", "ncRNA", 
   "novel antisense - filtered", "novel isolated - filtered", "unexpressed isolated")

##
## PIE: 4 classes
##
if("pie" %in% what){

  for(kindOfPie in 1:2) {
  
  if(interact) {
    x11(width=7*length(rnaTypes), height=4.8)
  } else {
    ## pdf("tableSegments-pie.pdf", width=7*length(rnaTypes), height=4.8)
    pdf(paste("tableSegments-pie-", kindOfPie, ".pdf", sep=""), width=7*length(rnaTypes), height=4.8)
  }

  par(mfrow=c(1, length(rnaTypes)))
  counts = NULL
  for(rt in rnaTypes) {
    s  = cs[[rt]]
    category = s[, "category"]
    overlap  = s[, "overlap"]

    px = switch(kindOfPie,
      c("overlap >=50%"   = sum(overlap  %in% c(">=50%, <100%", "100%")),
        "overlap < 50%"   = sum(overlap  %in% c("<50%")),
        "novel isolated - filtered"  = sum(category %in% c("novel isolated - filtered")),
        "novel isolated - excluded"  = sum(category %in% c("novel isolated - excluded")),
        "novel antisense - filtered" = sum(category %in% c("novel antisense - filtered")),
        "novel antisense - excluded" = sum(category %in% c("novel antisense - excluded"))),
      c("overlap >=50%"   = sum(overlap  %in% c(">=50%, <100%", "100%")),
        "overlap < 50%"   = sum(overlap  %in% c("<50%")),
        "novel isolated"  = sum(category %in% c("novel isolated - filtered", "novel isolated - excluded")),
        "novel antisense" = sum(category %in% c("novel antisense - filtered", "novel antisense - excluded"))),
      stop("Sapperlot"))
    
    stopifnot(all(names(px) %in% names(fillColors)))
    counts = cbind(counts, px)
    pie(px, radius=0.75, main=longNames[rt], 
        col = fillColors[names(px)],
        labels = paste(names(px), " (", px, ")", sep=""))

    if(kindOfPie==1){
      cat(">>>>>>>>>", rt, "<<<<<<<<<<<<\n")
      tab = table(category, overlap)
      tab = tab[rowSums(tab)!=0, ]
      print(tab)
      cat("\n\n")
    }
    
  }
  if(!interact)
    dev.off()

  cat(c("With", "Without")[kindOfPie], "Filtering:\n----------------------------\n")
  colnames(counts)=rnaTypes
  rownames(counts)=names(px)
  cat("\nSegment counts:\n")
  print(counts)
  cat("\n\n")
} ## for kindofPie
} ## if


##
## WRITE THE SEGMENT TABLE
##
if("wst" %in% what){
  selectedCategories = c(
     "annotated ORF", "ncRNA", "other annotation",
      "novel isolated - filtered", 
      "novel antisense - filtered")
  for(rt in rnaTypes) {
    s   = cs[[rt]]
    stopifnot(all(selectedCategories %in% levels(s[,"category"])))
    sel = (s[,"category"] %in% selectedCategories)
    fn  = file.path(indir[rt], "viz", "index.html")
    cat("Writing", fn, "\n")
    writeSegmentTable(s[sel, ], fn=fn, sortBy="category",
      title=paste(rt, " (", longNames[rt], ")", sep=""))
  }
  cat("\n")
}

##
## LENGTH DISTRIBUTIONS
##
maxlen=5000
if("length" %in% what){
  selectedCategories = c(
     "annotated ORF", "ncRNA", 
      "novel isolated - filtered", 
      "novel antisense - filtered")

  stopifnot(all(selectedCategories %in% names(fillColors)))
  if(!interact)
    pdf("tableSegments-lengths.pdf", width=14, height=length(rnaTypes)*3)
  par(mfrow=c(length(rnaTypes),5))
  br = seq(0, maxlen, by=200)
  for(rt in rnaTypes) {
    s   = cs[[rt]]
    stopifnot(all(selectedCategories %in% levels(s[,"category"])))
    for(lev in selectedCategories) {
      len = s[s[,"category"] == lev, "length"]
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
  selectedCategories = c("annotated ORF", "novel isolated - filtered")
  for(rt in rnaTypes) {
    s   = cs[[rt]]
    stopifnot(all(selectedCategories %in% levels(s[,"category"])))
    ylim = quantile(s[,"level"][s[,"category"] %in% selectedCategories],
      probs=c(0.01, 0.99), na.rm=TRUE)
    for(lev in selectedCategories) {
      len = s[s[,"category"] == lev, "length"]
      len[len>maxlen] = maxlen
      exl = s[s[,"category"] == lev, "level"]
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
    ## fas = br[[3]] * br[[4]] / s$length[br[[1]]]
    fas = (br[[4]] / s$length[br[[1]]] > 0.5) * 100 
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

  cat("Fraction of alignable sequence (percent):\n",
      "=========================================\n",
      "Here, 'number' is the number of segments.\n", sep="")

  theSplit = vector(mode="list", length=length(rnaTypes))
  names(theSplit) = rnaTypes
    
  selectedCategories = c("annotated ORF", "ncRNA",  "novel antisense - filtered", "novel isolated - filtered",
    "unexpressed isolated")

  for(rt in rnaTypes) {
    s   = cs[[rt]]

    catg = s[,"category"]
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
  cat("\n\nGrouping of conservation scores by expression, using ", nrlevs, "\n", 
      "quantile groups of equal size.\n",
      "=====================================================================\n", 
      sep="")
  
  if(interact) {
    x11(width=6*length(rnaTypes), height=8)
  } else {
    pdf(paste("tableSegments-consex.pdf", sep=""),
        width=4*length(rnaTypes), height=8)
  }
  par(mfcol=c(1, length(rnaTypes)))
  stopifnot(all(selectedCategories %in% names(lineColors)))

  for(i in seq(along=rnaTypes)) {
    rt=rnaTypes[i]
    cat("\n", rt, "\n------------\n", sep="")
    
    sp = theSplit[[rt]]
    s  = cs[[rt]]
        
    fraction = matrix(NA, nrow=nrlevs, ncol=length(sp))
    colnames(fraction) = names(sp)
  
    pchs = c(15, 5, 16, 17, 19) # 20
    stopifnot(length(pchs)==length(sp))
    
    for(j in seq(along=sp)) {
      wh   = sp[[j]]
      v    = s$level[wh]
      br   = quantile(v, (0:nrlevs)/nrlevs, na.rm=TRUE)
      spwh = split(wh, cut(v, breaks=br))
      hit  = calchit(spwh, blastres[[rt]], s)
      fraction[, names(sp)[j]] = hit$fraction
      cat("\nby expression for:", names(sp)[j], "\n")
      print(round(hit,1))
    }
    #j = which(colnames(fraction)=="isolated and unexpressed")
    #stopifnot(length(j)==1)
    #baseline = mean(fraction[, j])
    #fraction = fraction[, -j]
    
    matplot(fraction, xaxt="n", type="b", main=longNames[rt],
            lty=1, lwd=2, pch=pchs, col=lineColors[colnames(fraction)],
            ylab="average identity (percent) ", xlab="transcript level",
            ylim=c(0, 110))
    #abline(h=baseline, col=lineColors[["isolated and unexpressed"]], lwd=2, type="l")
    cutNames = paste("<=", round((1:nrlevs)/nrlevs*100), "%", sep="")
    axis(side=1, at = 1:nrlevs, labels = cutNames)
    if(i==1)
      legend(x=1, y=112, legend=colnames(fraction),
             lty=1, lwd=2, pch=pchs, col=lineColors[colnames(fraction)])
  }
  if(!interact)
    dev.off()
 
}

if(!interact)
  sink()
