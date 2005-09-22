##
## How much of the genome is transcribed
## 1. By ORFs
## 2. By Basepairs --> (Fig.4)
##
## As background, we use in both cases the background estimated from the segmentation (in calcThreshold).
## For the ORF detection, we detect in poly-A and total separately, using FDR threshold 0.1%,
## then take the union of the result.
##

library("tilingArray")
library("multtest")
source("setScriptsDir.R")

rnaTypes  = c("seg-polyA-050909", "seg-tot-050909")
source(scriptsDir("readSegments.R"))

interact=!TRUE
if(!interact) {
  sink("expressedFeatures.txt")
  cat("Made on", date(), "\n\n")
}

source(scriptsDir("calcThreshold.R"))


## Load the raw data
for(rt in rnaTypes) {
  e = get(rt)
  if(!("xn" %in% ls(e))) {
    fn = file.path(rt, "xn.rda")
    cat("Loading", fn, "\n")
    load(fn, envir=e)
    x = exprs(get("xn", e))
    x = rowMeans(x) - get("theThreshold", e)
    assign("x", x, envir=e)
    rm(x)
  }
}


sel = gff[, "feature"]=="gene"
allncRNA = c("ncRNA","snoRNA","snRNA","tRNA","rRNA")

if(!"ID"%in%colnames(gff))
  gff$ID=getAttributeField(gff[, "attributes"], "ID")

sgff = rbind(
  cbind(category="verified gene",        gff[ sel & gff[, "orf_classification"]=="Verified", ]),
  cbind(category="uncharacterized gene", gff[ sel & gff[, "orf_classification"]=="Uncharacterized", ]),
  cbind(category="dubious gene",         gff[ sel & gff[, "orf_classification"]=="Dubious", ]),
  cbind(category="ncRNA(all)",           gff[ gff[, "feature"] %in% allncRNA & !is.na(gff[, "ID"]), ]))
colnames(sgff)[1]="category"

cat("Features in GFF table:\n")
print(table(sgff$category))


stopifnot(!any(is.na(sgff[, "Name"])), !any(duplicated(sgff[, "Name"])))

## for each interesting feature (gene, RNA), build a list of probes:
if(!exists("indProbe")) {
  indProbe = vector(mode="list", length=nrow(sgff))
  names(indProbe) = sgff[, "Name"]
  
  cat("Calculating indProbe (length=", length(indProbe), "): ", sep="")

  nm.ind = paste(gff$chr, gff$strand, "index", sep=".")
  nm.sta = paste(gff$chr, gff$strand, "start", sep=".")
  nm.end = paste(gff$chr, gff$strand, "end",   sep=".")
  nm.uni = paste(gff$chr, gff$strand, "unique",   sep=".")
  
  sp = split(1:nrow(gff), gff[, "Name"])
  for(i in seq(along=indProbe)) {
    if(i%%250==0) cat(i, "")
    ## this combination of for/if is really slow
    whf = sp[[names(indProbe)[i]]]
    fn  = gff[whf, "feature"]
    ## genes (may contain introns)
    if(sum(fn=="gene")==1 && all(fn%in%c("gene", "CDS", "intron", "region"))) {
      fset = which(fn=="CDS")
      ## tRNAs (may contain introns)  
    } else if (sum(fn=="tRNA")==1 && all(fn%in%c("tRNA","ncRNA","intron"))) {
      fset = which(fn=="ncRNA")
      ## rRNAs (may contain introns)  
    } else if (sum(fn=="rRNA")==1 && all(fn%in%c("rRNA","ncRNA","intron", "nc_primary_transcript"))) {
      fset = which(fn%in%c("ncRNA","nc_primary_transcript"))
      ## snoRNAs (may contain introns)  
    } else if (sum(fn=="snoRNA")==1 && all(fn%in%c("snoRNA","ncRNA", "intron"))) {
      fset = which(fn=="ncRNA")
      ## snRNAs   
    } else if (sum(fn=="snRNA")==1 && all(fn%in%c("snRNA","ncRNA"))) {
      fset = which(fn=="ncRNA")
    } else if (all(fn=="ncRNA")) {
      fset = which(fn=="ncRNA")
    } else {
      ##cat("Dropping:\n")
      ##print(gff[whf, c(1, 3:5, 6)])
      ##cat(gsub("%20", " ", unique(getAttributeField(gff[whf, "attributes"], "Note"))), "\n\n")
      fset = integer(0)
    }
    ## gene with one or more CDSs
    res = integer(0)
    for(w in whf[fset])
      res = c(res, get(nm.ind[w], probeAnno)[
        (get(nm.uni[w], probeAnno) == 0) &
        (get(nm.sta[w], probeAnno) >= gff$start[w])&
        (get(nm.end[w], probeAnno) <= gff$end[w]) ])
    indProbe[[i]] = sort(unique(res))
  }
  cat("\n")
}

## select features that have >= 7 probes
hasEnoughProbes = which(listLen(indProbe)>=7)
cat(length(hasEnoughProbes),"of",length(indProbe),"potential transcripts are matched by >=7 unique probes.\n\n")

## do the testing
stopifnot(all(rnaTypes==c("seg-polyA-050909", "seg-tot-050909")))
res   = matrix(NA, nrow=length(levels(sgff[,"category"])), ncol=11)
isExp = matrix(as.logical(NA), nrow=length(hasEnoughProbes), ncol=length(rnaTypes)+1)

rownames(res) = levels(sgff[,"category"])
colnames(res) = c("n1: in genome", "n2: with probes", 
    "detected in poly-A RNA", "% of n1", "% of n2",
    "detected in total RNA",  "% of n1", "% of n2",
    "detected in either", "% of n1",  "% of n2")

res[names(n1), 1] = n1 = table(sgff[, "category"])
res[names(n2), 2] = n2 = table(sgff[hasEnoughProbes, "category"])


for(irt in 1:3) {

  if(irt<=2){
    rt = rnaTypes[irt]
    x = get("x", get(rt))

    ## The reason behind this is a little bit more complicated than it seems:
    ## 'x' is the average of 3 arrays (in poly-A case) resp. 2 arrays (total RNA case).
    ## Hence we seem to loose power by doing the sign test on the averages rather than
    ## on the individual values (which across arrays are independent).
    ## On the other hand side, the probes overlap: their length is 25, their start sites
    ## typically 8 bases apart. By ignoring the correlation between neighbouring probes,
    ## the sign test will reject to easily. The assumption behind this is that two effects
    ## roughly cancel; moreover, that for most probe sets the real effect is so overwhelming
    ## that these comparatively small effects do not matter.
    doTest = function(ind) {
      pv = sapply(ind, function(jj)
        binom.test(sum(x[jj]>0), length(jj), alternative="greater")$p.value)
      
      bh = mt.rawp2adjp(pv, proc="BY")
      stopifnot(all(bh$adjp[, 1] == pv[bh$index], na.rm=TRUE))
      adjp = numeric(nrow(bh$adjp))
      adjp[bh$index] = bh$adjp[,2]
      (adjp < FDRthresh)
    }

    isExp[, irt] = doTest(indProbe[hasEnoughProbes])
  } else {
    isExp[, irt] = (isExp[, 1] | isExp[, 2])
  }
    
  nrDet = table(sgff[hasEnoughProbes, "category"], isExp[, irt])[, "TRUE"]

  stopifnot(identical(names(nrDet), names(n1)), identical(names(nrDet), names(n2)))
  stopifnot(setequal(names(nrDet), rownames(res)))
  
  res[names(nrDet), irt*3]   = nrDet
  res[names(nrDet), irt*3+1] = round(100*nrDet/n1, 1)
  res[names(nrDet), irt*3+2] = round(100*nrDet/n2, 1)    
}
  
print(res)
cat("\n\n\n\n")

##
## What fraction of probes in the genome are transcribed
##

data(yeastFeatures)
transcribedFeatures = rownames(yeastFeatures)[yeastFeatures$isTranscribed]

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
  s = get("segScore", get(rt))
  lev = s[, "level"]
  res = lapply(1:nrChr, function(chr) {
    res  = rep(-Inf, chrlen[chr])
    selt = which(s[, "chr"]==chr & !is.na(s[, "level"]) & s[,"frac.dup"]<maxDuplicated)
    for(i in selt) {
      rg = s$start[i]:s$end[i] 
      res[rg] = pmax(res[rg], lev[i])
    }
    res
  })
  segLev[[rt]]  = unlist(res)
  isTrans[[rt]] = (segLev[[rt]] >= 0)
  stopifnot(length(segLev[[rt]]) == length(isAnno))
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
  pdf("expressedFeatures.pdf", height=3, width=4)
par(mfrow=c(1,1))
myHist = function(x) {
  xmax = quantile(x, 0.9999, na.rm=TRUE)
  xmin = min(x[is.finite(x)])
  x[x>xmax] = xmax
  by = 0.1
  breaks = c(rev(seq(0, xmin-by, by=-by)), seq(by, xmax+by, by=by))
  theCol = brewer.pal(4, "Paired")[3]
  ##hist(x, breaks=breaks, col=cols[1], main="", yaxt="n", ylab="", xlab="level")
  sf = showDens(z=list(x=x), breaks=breaks, col=theCol, main="",  xlab="expression level")
  
  axis(side=2, at=(0:2)*2e5*sf[1], labels=c("0", "200000", "400000"), las=1)
  abline(v=0, col="black", lwd=3)
}
myHist(segLev[["both"]])

if(!interact)
  dev.off()  


##
##  GO analysis of unexpressed genes
##
cat("\n\nGO-Analysis of the untranscribed verified genes:\n\n")

source(scriptsDir("GOHyperG.R"))
source(scriptsDir("writeSegmentTable.R"))

unexpressedGenes = names(indProbe)[hasEnoughProbes][!isExp[,3] & sgff[hasEnoughProbes, "category"]=="verified gene"]
GOHyperG(unexpressedGenes)

if(!interact)
  sink()





