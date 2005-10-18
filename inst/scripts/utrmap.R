## 5 plots:
## a) distribution of 5' UTR lengths
## b) distribution of 3' UTR lengths
## c) scatterplot 5' vs 3' UTR lengths
## d) scatterplot mean level vs 5' UTR length
## e) scatterplot mean level vs 3' UTR length
## The poly-A version is for the paper,
## the total RNA version is for the supplement.

library("tilingArray")
library("davidTiling")
library("geneplotter")
source("setScriptsDir.R")

interact=(!TRUE)
options(error=recover, warn=0)
graphics.off()

if(!interact) {
  sink("utrmap.txt")
}

rnaTypes  =  c("seg-polyA-050909", "seg-tot-050909")
source(scriptsDir("readSegments.R"))
source(scriptsDir("calcThreshold.R"))
source(scriptsDir("categorizeSegments.R")) 
source(scriptsDir("writeSegmentTable.R"))

what = c("stat", "wst", "explen", "polyAvstot", "go")[-2]

##
## CATEGORIZE
##
if(!exists("cs")) {

  utr = vector(mode="list", length=length(rnaTypes)+1)
  names(utr) = c(rnaTypes, "combined")
  cs  = vector(mode="list", length=length(rnaTypes))
  names(cs) = rnaTypes

  for(rt in rnaTypes) {
    cat("\n--------", rt, "---------\n")
    s = categorizeSegmentsUTRmap(get(rt))
    s = s[!is.na(s[,"goodUTR"]), ]
    z = as.matrix(s[, c("utr5", "utr3")])
    rownames(z) = rownames(s) = s[, "featureInSegment"]
    colnames(z) = c("5' UTR", "3' UTR")
    utr[[rt]] = z
    cs[[rt]] =s
  }
  
  utr[["combined"]] = rbind(utr[[1]], utr[[2]][ !(rownames(utr[[2]])%in%rownames(utr[[1]])), ])
  save(utr, file=paste("utr-", date(), ".rda", sep=""))
  rm(list=c("s", "z"))
  
} else {
  cat("\n**************************************************\n",
        "*      NOT REDOING categorizeSegments            *\n",
        "**************************************************\n", sep="")
}

graphics.off()
cols = brewer.pal(12, "Paired")

scatterWithHist = function(x, breaks, barcols, xlab, ylab, ...) {
  stopifnot(is.matrix(x), ncol(x)==2)
  xmax = breaks[length(breaks)]
  xmid = breaks[length(breaks)/2]
  x[x>xmax]=NA
  xhist = hist(x[,1], breaks=breaks, plot=FALSE)
  yhist = hist(x[,2], breaks=breaks, plot=FALSE)
  topx  = max(xhist$counts)
  topy  = max(yhist$counts)
  top   = max(topx, topy)
  xrange = yrange = breaks[c(1, length(breaks))]
  nf = layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
 
  par(mar=c(3,3,1,1))
  plot(x, xlim=xrange, ylim=yrange, xlab="", ylab="", ...)
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0, col=barcols[1])
  text(length(xhist$counts)/2, topx, adj=c(0.5, 1), labels=xlab, cex=1.6)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, col=barcols[2], horiz=TRUE)
  text(topy, length(yhist$counts)/2, adj=c(0.5, 1), labels=ylab, cex=1.6, srt=270)
} 

##
## distribution summaries, length histograms, scatterplot 3' vs 5' length
##
if("stat" %in% what){
  for(rt in rnaTypes) {
    ul = utr[[rt]]
    cat("\n-------", rt, "-------\n")
    cat("Length distribution summary of", nrow(ul), "5'-UTRs:\n")
    print(summary(ul[, "5' UTR"]))
    cat("Length distribution summary of", nrow(ul), "3'-UTRs:\n")
    print(summary(ul[, "3' UTR"]))

    if(!interact) {
      pdf(file=paste("utrmap-", rt, ".pdf", sep=""), width = 5.5, height = 5.5)
    } else {
      x11(width = 7, height = 7)
    }
    maxlen = 700
    br = seq(0, maxlen, by=20)
    scatterWithHist(ul,
         xlab=paste("Length of", colnames(ul)[1]),
         ylab=paste("Length of", colnames(ul)[2]),
         breaks = br, pch=20, barcols=cols[c(1,3)])
    
    if(!interact)
      dev.off()
  }

  ## common:
  comUTR = intersect(rownames(utr[[1]]), rownames(utr[[2]]))
  allUTR = union(rownames(utr[[1]]), rownames(utr[[2]]))
  cat("\n", length(comUTR)," in both ", paste(rnaTypes, collapse=" and "), ", ", 
      length(allUTR), " altogether.\n", sep="")
}


##
## WRITE THE SEGMENT TABLE
##
if("wst" %in% what){
  for(rt in rnaTypes) {
    fn = file.path(indir[rt], "viz", "utrmap")
    nr = nrow(cs[[rt]])
    if(interact)
      cat("Writing", nr, "UTRs to", fn, "\n")
    writeSegmentTable(cs[[rt]], title=paste(nr, "UTR maps from", longNames[rt]), fn=fn,
                      sortBy = "goodUTR", sortDecreasing=TRUE, interact=interact)
  }
}

##
## expression vs length
##
if("explen" %in% what){

  investigateExpressionVersusLength = function(lev, len, main) {
    theCut = cut(len, breaks=quantile(len, probs=c(0, 0.95, 1)))
    e = tapply(lev, theCut, ecdf)
    theCol = cols[1]
    plot(e[[1]], pch=".", xlab="level", main=main)
    for(i in 2:length(e)) {
      theCol = cols[i*2]
      lines(e[[i]], col.hor=theCol, col.points=theCol, col.vert=theCol)
    }
  }

  investigateLengthVersusLength = function(l1, l2, ...) {
    cc = cor(l1, l2, method="kendall")
    plot(l1+1, l2+1, log="xy", pch=".", main=paste("length, cor=", signif(cc, 3)), ...)
  }
  
  if(!interact) {
    pdf(file=paste("utrmap-expression-vs-length.pdf", sep=""), width = 10.5, height = 14)
  } else {
    x11(width = 10.5, height = 14)
  }
  par(mfrow = c(4, 3))

  for(rt in rnaTypes) {
    s = cs[[rt]]
    ## get the CDS length
    mt = match(s[,"featureInSegment"], gff[,"Name"])
    stopifnot(!any(is.na(mt)))
    cdslen = gff[mt, "end"]-gff[mt, "start"]

    investigateExpressionVersusLength(s[,"level"], s[,"utr3"], paste(rt, ": length of 3' UTR", sep=""))
    investigateExpressionVersusLength(s[,"level"], s[,"utr5"], paste(rt, ": length of 5' UTR", sep=""))
    investigateExpressionVersusLength(s[,"level"], cdslen, paste(rt, ": length of CDS", sep=""))
    investigateLengthVersusLength(s[,"utr3"], s[,"utr5"], xlab="3' UTR", ylab="5' UTR")
    investigateLengthVersusLength(cdslen, s[,"utr3"], xlab="CDS", ylab="3' UTR")
    investigateLengthVersusLength(cdslen, s[,"utr5"], xlab="CDS", ylab="5' UTR")
  }
  
  if(!interact)
    dev.off()
}

##
## difference between total and poly-A
## 
if("polyAvstot" %in% what){
  if(interact) {
    x11(width=6.6, height=10)
  } else {
    pdf(file=paste("utrmap-scatter.pdf", sep=""), width=6.6, height=10)
  }
  
  par(mfrow=c(3,2))
  for(i in 1:2){
    px = utr[[1]][comUTR,i]
    py = utr[[2]][comUTR,i]
    axlim = c(0, quantile(c(px, py), 0.8))
    plot(px, py,
         main = paste("length of ", colnames(utr[[1]])[i], " (", length(comUTR), " common)", sep=""),
         xlab = longNames[rnaTypes[1]], ylab=longNames[rnaTypes[2]],
         xlim = axlim, ylim = axlim, pch=20)
    abline(a=0, b=1, col="red")
  }

  vec = c("5' UTR", "3' UTR")
  d = utr[[1]][comUTR, vec] - utr[[2]][comUTR, vec]
  colnames(d)=vec
  
  ex1 = cs[[1]][comUTR, "level"]
  ex2 = cs[[2]][comUTR, "level"]

  ex = cbind(difference=ex1-ex2, average=(ex1+ex2)/2)

  col2 = "lightblue"; col3="pink"
  for(j in 1:ncol(ex)) {
    for(i in 1:ncol(d)){
      ec1 = ecdf(ex[d[,i]==0, j])
      ec2 = ecdf(ex[d[,i]>0, j])
      ec3 = ecdf(ex[d[,i]<0, j])
      plot(ec1, pch=".", xlab=colnames(ex)[j], main=colnames(d)[i])
      lines(ec2, pch=".", col.points=col2, col.hor=col2, col.vert=col2)
      lines(ec3, pch=".", col.points=col3, col.hor=col3, col.vert=col3)
    }
  }
  if(!interact)
    dev.off()
}

##
## GO analysis of UTR lengths
##
if("go" %in% what){

  ## Use poly-A data only
  myUTR = utr[["seg-polyA-050909"]]  
  goCat = getAllGO(rownames(myUTR), gff)
  
  allGO = unique(unlist(goCat))
  gm    = matrix(FALSE, nrow=length(allGO), ncol=length(goCat))
  rownames(gm) = allGO
  colnames(gm) = names(goCat)
  
  for(i in seq(along=goCat)) {
    gm[ goCat[[i]], i] = TRUE
  }
  w.all = which(rownames(gm)=="all")
  stopifnot(length(w.all)==1, all(gm[w.all, ]))
  gm = gm[-w.all,]
  
  wtfun = function() {
    function(z) {
      sz = sum(z)
      if(sz>=5&&(length(z)-sz)>=5) {
        w5 = wilcox.test(myUTR[, "5' UTR"] ~ z)$p.value
        w3 = wilcox.test(myUTR[, "3' UTR"] ~ z)$p.value
      } else {
        w5 = w3 = as.numeric(NA)
      }
      m5 = median(myUTR[z, "5' UTR"])
      m3 = median(myUTR[z, "3' UTR"])
      c(w5, m5, w3, m3, sz)
    }
  }
  
  wt = wtfun()
  pGO = t(apply(gm, 1, wt))
  colnames(pGO) = c("p 5' UTR", "median 5' UTR", "p 3' UTR", "median 3' UTR", "nrGenes")

  medAll = apply(myUTR, 2, median)
  
  ## print the GO TERMS
  GOterms = mget(rownames(gm), GOTERM)

  for(j in 1:ncol(myUTR)) {
    pname = paste("p",      colnames(myUTR)[j])
    mname = paste("median", colnames(myUTR)[j])
    cat("--------------------------------------------------\n",
        colnames(myUTR)[j], "\n",
        "--------------------------------------------------\n", sep="")
    sel = order(pGO[, pname])[1:40]
    for(s in sel) {
      g   = rownames(gm)[s]
      cat(g, " p=", format.pval(pGO[s, pname]), "  median=",
          pGO[s, mname], " (versus ", medAll[j], ")\n", sep="")
      print(GOterms[[s]])
      nr = pGO[s, "nrGenes"]
      cat(nr, " genes", sep="")
      if(nr<=20)
        cat(":", replaceSystematicByCommonName(names(which(gm[s,]))))
      cat("\n\n")
    }
  }

  ## Write tab-delimited table
  ## GO-ID, Term, Median 3' UTR, Median 5' UTR, No. Genes
  pp  = pmin(pGO[, "p 5' UTR"], pGO[, "p 3' UTR"])
  ord = order(pp)
  ord = ord[ which(pp[ord]<= 2e-3) ]
  sGOterms = GOterms[ord]
  
  outtab = cbind(data.frame(
    GOID = sapply(sGOterms, GOID),
    Term = sapply(sGOterms, Term), 
    Ontology = sapply(sGOterms, Ontology)), 
    pGO[ord, ])

  write.table(outtab, file="utrmap-GOcategories.txt", sep="\t", row.names=FALSE)

  
  
  ## Look at 3' UTR lengths of Cellular Components
  cat("-----------------------------------------------------------------\n",
      "Cellular Component Categories with large median length of 3' UTRs\n",
      "-----------------------------------------------------------------\n\n\n", sep="")

  cellularComponents = names(GOterms)[ sapply(GOterms, Ontology) == "CC" ]

  threePrimeUTRLengths = vector(mode="list", length=length(cellularComponents))
  names(threePrimeUTRLengths) = cellularComponents
  for(i in seq(along=cellularComponents)) {
    genes = colnames(gm)[ gm[cellularComponents[i], ] ]
    tmp = myUTR[ genes, "3' UTR"]
    names(tmp) = genes
    threePrimeUTRLengths[[i]] = tmp
  }

  
  ord = order(sapply(threePrimeUTRLengths, median), decreasing=TRUE)
  for(j in ord[1:25]) {
    print(get(names(threePrimeUTRLengths)[j], GOTERM))
    print(threePrimeUTRLengths[[j]])
    cat("\n\n")
  }  
}

if(!interact) {
  sink()
}

