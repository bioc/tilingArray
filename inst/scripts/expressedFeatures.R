##
## How many known features do we find exprssed?
##
library("tilingArray")
library("multtest")
source("scripts/readSegments.R") 

interact=!TRUE
if(!interact) {
  sink("expressedFeatures.txt")
  cat("Made on", date(), "\n\n")
}

source("scripts/calcThreshold.R") 

if(!exists("intergenic")){
  probe = probeAnno$probeReverse
  intergenic = which(probeAnno$probeReverse$no_feature=="no" & probeAnno$probeDirect$no_feature=="no")
}

for(rt in rnaTypes) {
  e = get(rt)
  if(!("xn" %in% ls(e))) {
    fn = file.path(rt, "xn.rda")
    cat("Loading", fn, "\n")
    load(fn, envir=e)
    x = exprs(get("xn", e))
    for(j in 1:ncol(x))
      x[,j] = x[, j] - median(x[intergenic, j])
    assign("x", x, envir=e)
    rm(x)
  }
}

testExp = function(index) {
  pv = sapply(index, function(jj) {
    rv = 1
    if(length(jj)>0)
      rv = binom.test(sum(x[jj, ]>0), ncol(x)*length(jj), alternative="greater")$p.value
    rv
  })
    
  bh = mt.rawp2adjp(pv, proc="BY")
  stopifnot(all(bh$adjp[, 1] == pv[bh$index], na.rm=TRUE))
  adjp = numeric(nrow(bh$adjp))
  adjp[bh$index] = bh$adjp[,2]
  (adjp < FDRthresh)
}

sel = gff[, "feature"]=="gene"
allncRNA = c("ncRNA","snoRNA","snRNA","tRNA","rRNA")
categGff = list(
  "verified gene"        = gff[ sel & gff[, "orf_classification"]=="Verified", ],
  "uncharacterized gene" = gff[ sel & gff[, "orf_classification"]=="Uncharacterized", ],
  "dubious gene"         = gff[ sel & gff[, "orf_classification"]=="Dubious", ],
  "ncRNA(all)"           = gff[ gff[, "feature"] %in% allncRNA , ])
stopifnot(all(sapply(categGff, nrow)>=1))


res = matrix(NA, nrow=length(categGff), ncol=8)
stopifnot(all(rnaTypes==c("seg-polyA-050525", "seg-tot-050525")))
rownames(res) = names(categGff)
colnames(res) = c("n1: in genome", "n2: with probes", 
    "det. poly-A", "poly-A: % of n1", "poly-A: % of n2",
    "det. total",  "total: % of n1",  "total: % of n2")

for(irt in seq(along=rnaTypes)) {
  rt = rnaTypes[irt]
  cat("\n=====", rt, "=====\n")
  x = get(rt)$"x"
    
  for(icg in seq(along=categGff)) {
    sgff = categGff[[icg]]

    ## build a list of probes
    unames = sort(unique(sgff$Name))
    stopifnot(!any(unames==""))
    index = vector(mode="list", length=length(unames))
    names(index) = unames

    ind = paste(sgff$chr, sgff$strand, "index", sep=".")
    sta = paste(sgff$chr, sgff$strand, "start", sep=".")
    end = paste(sgff$chr, sgff$strand, "end",   sep=".")
    uni = paste(sgff$chr, sgff$strand, "unique",   sep=".")

    for(i in 1:nrow(sgff)) {
      nm = sgff$Name[i]
      index[[nm]] = c(index[[nm]], get(ind[i], probeAnno)[
       get(uni[i], probeAnno) &
      (get(sta[i], probeAnno) >= sgff$start[i])&
      (get(end[i], probeAnno) <= sgff$end[i]) ])
    }

    ## How many features have >= 7 probes
    n1 = length(index)
    n2 = sum(listLen(index) >= 7)
    if(irt==1) {
      res[icg, 1:2] = c(n1,n2)
    } else {
      stopifnot(all(res[icg, 1:2] == c(n1,n2)))
    }
    
    ## Test
    isExp = testExp(index)
    nrDet = sum(isExp)
    res[icg, irt*3] = nrDet
    res[icg, irt*3+1] = round(100*nrDet/n1, 1)    
    res[icg, irt*3+2] = round(100*nrDet/n2, 1)    

    ## controls
    ctrlIndex = lapply(index, function(x) sample(intergenic, length(x)))
    ctrlExp = testExp(ctrlIndex)
    cat(sprintf("%25s: %4d, control: %4d\n", 
        names(categGff)[icg], as.integer(nrDet), as.integer(sum(ctrlExp))))
  }
  
}

print(res)

if(!interact)
  sink()
