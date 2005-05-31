##
## How many known features do we find exprssed?
##
library("tilingArray")
library("multtest")
source("scripts/readSegments.R") 
fdrThresh = 0.001

interact=!TRUE
if(!interact) {
  sink("expressedFeatures.txt")
  cat("Made on", date(), "\n\n")
}

cat("fdrTresh=", fdrThresh, "\n")
if(!exists("intergenic")){
  probe = probeAnno$probeReverse
  intergenic = which(probeAnno$probeReverse$no_feature=="no" & probeAnno$probeDirect$no_feature=="no")
}

if(!exists("allxn")) {
  allxn = vector(mode="list", length=length(rnaTypes))
  names(allxn) = rnaTypes
  for(rt in rnaTypes) {
    fn = file.path(rt, "xn.rda")
    cat(fn, "")
    load(fn)
    x = exprs(xn)
    for(j in 1:ncol(x))
      x[,j] = x[, j] - median(x[intergenic, j])
    allxn[[rt]] = x 
  }
  cat("\n")
}



testExp = function(index) {
  pv = sapply(index, function(ind)
    binom.test(sum(x[ind, ]>0), ncol(x)*length(ind), alternative="greater")$p.value
    )
    
  bh = mt.rawp2adjp(pv, proc="BY")
  stopifnot(all(bh$adjp[, 1] == pv[bh$index], na.rm=TRUE))
  adjp = numeric(nrow(bh$adjp))
  adjp[bh$index] = bh$adjp[,2]
  (adjp < fdrThresh)
}

data(transcribedFeatures)
feats = transcribedFeatures[transcribedFeatures!="gene"]

for(rt in rnaTypes) {
  cat("\n=====", rt, "=====\n")
  x = allxn[[rt]]
  
  for(ft in feats) {

    sgff = gff[ gff[, "feature"]==ft, ]
    stopifnot(nrow(sgff)>=1)

    ## build a list of probes
    unames = sort(unique(sgff$Name))
    stopifnot(!any(unames==""))
    index = vector(mode="list", length=length(unames))
    names(index) = unames

    ind = paste(sgff$chr, sgff$strand, "index", sep=".")
    sta = paste(sgff$chr, sgff$strand, "start", sep=".")
    end = paste(sgff$chr, sgff$strand, "end",   sep=".")
    for(i in 1:nrow(sgff)) {
      ## if(i%%100==0) cat(i, "")
      nm = sgff$Name[i]
      index[[nm]] = c(index[[nm]], get(ind[i], probeAnno)[
      (get(sta[i], probeAnno) >= sgff$start[i])&
      (get(end[i], probeAnno) <= sgff$end[i]) ])
    }
    
    isExp = testExp(index)

    ## controls
    ctrlIndex = lapply(index, function(x) sample(intergenic, length(x)))
    ctrlExp = testExp(ctrlIndex)
    
    cat(sprintf("%25s: %4d of %4d (%3.1f percent), control: %4d\n", 
                ft, as.integer(sum(isExp)), length(isExp), signif(100*mean(isExp)),
                as.integer(sum(ctrlExp))))
  }
  
}

if(!interact)
  sink()
