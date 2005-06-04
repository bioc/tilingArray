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

for(rt in rnaTypes) {
  cat("\n=====", rt, "=====\n")
  x = allxn[[rt]]
  
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
      ## if(i%%100==0) cat(i, "")
      nm = sgff$Name[i]
      index[[nm]] = c(index[[nm]], get(ind[i], probeAnno)[
       get(uni[i], probeAnno) &
      (get(sta[i], probeAnno) >= sgff$start[i])&
      (get(end[i], probeAnno) <= sgff$end[i]) ])
    }
    
    isExp = testExp(index)

    ## controls
    ctrlIndex = lapply(index, function(x) sample(intergenic, length(x)))
    ctrlExp = testExp(ctrlIndex)
    
    cat(sprintf("%25s: %4d of %4d (%3.1f percent), control: %4d\n", 
                names(categGff)[icg], as.integer(sum(isExp)), length(isExp), signif(100*mean(isExp)),
                as.integer(sum(ctrlExp))))
  }
  
}

if(!interact)
  sink()
