library("vsn")
library("genefilter")

if(!exists("a"))load("a.rda")
if(!exists("probeAnno"))load("probeAnno.rda")

rrr = 5
hybeSets = list(
  "polyA2" = c("05_04_27_2xpolyA_NAP3.cel.gz",
    "05_04_26_2xpolyA_NAP2.cel.gz",
    "05_04_20_2xpolyA_NAP_2to1.cel.gz"),
  "tot" = c("050409_totcDNA_14ug_no52.cel.gz",
    "030505_totcDNA_15ug_affy.cel.gz"),
  "tot2" = c("050331_totcDNA_15ug_S96.cel.gz",
    "050411_totcDNA_20ug_affy.cel.gz",
    "050415_totcDNA_20ug_Affy11.cel.gz"),
  "dir" = c("050621_dirPolyARNA_10ug_2-3.cel.gz",
    "050621_dirPolyARNA_10ug_2-3_4x.cel.gz"),
  "odT" = c("041112_S96_polyA-dT-cDNA1_16H_45C.cel.gz"))[rrr]

outdir = c("polyA2" = "seg-polyA-050521",
           "tot"    = "seg-tot-050521",
           "tot2"   = "seg-tot2-050521",
           "dir"    = "seg-dir-050721",
           "odT"    = "seg-odT-050801")[rrr]

normMethod = c(rep("vsn", 4), "shiftlog")[rrr]
names(normMethod) = names(outdir)

##
## check
##
stopifnot(length(outdir)==length(hybeSets),
          names(hybeSets)==names(outdir),
          length(hybeSets)==length(normMethod),
          names(hybeSets)==names(normMethod))
          
##
## DNA hybes
## 
jref = which(a$NucleicAcid == "DNA")
stopifnot(length(jref)==3)

if(!exists("allPM")) {
  ## 1. select PM probes
  chrstr = paste(rep(1:17, each=2), 
    rep(c("+", "-"), 17), sep=".")
  allPM = unique(unlist(lapply(chrstr, function(chr)
    get(paste(chr, "index", sep="."), probeAnno))))
  
  ## 2. select intergenic probes
  allPM = allPM[ probeAnno$probeReverse$no_feature[allPM]=="no" &
    probeAnno$probeDirect$no_feature[allPM]=="no" ]
  
  cat("Selected", length(allPM), "intergenic PM probes.\n")
}

if(!exists("refSigPM")) {
  ## DNA-normalization factor:
  refSig    = rowMeans(log(exprs(a)[, jref, drop=FALSE], 2))
  refSigPM  = refSig[allPM]
  
  nrStrata  = 20
  quants    = quantile(refSig, probs=seq(0, 1, length=nrStrata+1))
  quants[1] = quants[1]-1
  strata    = cut(refSigPM, quants)
}

for(hs in seq(along=hybeSets)) {
  cat(names(hybeSets)[hs], "\n")
  fn = hybeSets[[hs]]
  
  x.bg  = tapply(refSigPM, strata, median)
  y.bg  = matrix(NA, nrow=length(x.bg), ncol=length(fn))
  bgfun = vector(mode="list", length=length(fn))
  colnames(y.bg) = names(bgfun) = fn
  for(f in fn) {
    y.bg[,f]   = tapply(log(exprs(a)[allPM, f], 2), strata, shorth)
    bgfun[[f]] = approxfun(x.bg, y.bg[,f], rule=2)
  }    
  matplot(x.bg, y.bg, pch=15+seq(along=fn))
  px = seq(5, 13, by=0.1)
  for(bf in bgfun)
    lines(px, bf(px))
  
  yn = matrix(NA, nrow=nrow(exprs(a)), ncol=length(fn))
  colnames(yn) = fn
  for(f in fn)
    yn[, f] = (exprs(a)[, f] - 2^bgfun[[f]](refSig) ) / 2^refSig

  ## vsn with strata? Not sure whether this is useful
  ## vsr = vsn(yn, lts.quantile=0.95, strata=as.integer(strata), subsample=1e5)

  switch(normMethod[hs],
       vsn = {
         ## vsn without strata
         xn = vsn(yn, lts.quantile=0.95, subsample=2e5)
         exprs(xn)    = exprs(xn)/log(2)
         phenoData(xn) = phenoData(a)[fn, ]
       },
       shiftlog = {
         offset = 0.5
         stopifnot(!any(yn+offset <= 0))
         xn = new("exprSet", exprs = log(yn+offset, 2), phenoData=phenoData(a)[fn, ])
       },
       stop("Zapperlot")
  )

  save(xn, refSig, file=file.path(outdir[hs], "xn.rda"), compress=TRUE)
}
