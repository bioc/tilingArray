library("Biobase")
library("RColorBrewer")
source(scriptsDir("readSegments.R"))
options(error=recover)

interact=TRUE
makeFig=!TRUE

if(!exists("a")) {
  cat("Loading a.rda\n")
  load("a.rda")
}

for(rt in rnaTypes) {
  if(!("xn" %in% ls(get(rt)))) {
    fn = file.path(rt, "xn.rda")
    cat("Loading", fn, "\n")
    load(fn, envir=get(rt))
  }
}

refSig = get("refSig", get(rnaTypes[1]))
for(j in seq(along=rnaTypes)[-1])
  stopifnot(identical(refSig, get(rnaTypes[j])$refSig))


if(!exists("refSigThresh")) {
  chrstr = paste(rep(1:17, each=2),
    rep(c("+", "-"), 17), sep=".")
  allPM = unique(unlist(lapply(chrstr, function(chr)
    get(paste(chr, "index", sep="."), probeAnno))))
  refSigThresh = quantile(refSig[allPM], probs=0.05)
}
cat("refSigThresh=", signif(refSigThresh, 3), "\n")

jref = which(a$NucleicAcid == "DNA")
stopifnot(length(jref)==3)

##
## the five different normalization methods
##
calcY = function(isel, fn, rt) {
  y1 = rowMeans(log(exprs(a)[isel, jref], 2))
  y2 = rowMeans(log(exprs(a)[isel, fn], 2))
  y3 = rowMeans(exprs(get(rt)$xn)[isel,  ])
  y4 = y3
  y4[refSig[isel]<refSigThresh] = NA
  y = cbind(y1, y2, y2-y1, y3, y4)
  colnames(y) = c("DNA", "unnormalized", "normalization method 1",
            "normalization method 2","normalization method 3")
  return(y)
}

calcdd = function(x, k=4)
 abs(x[-(1:k)] - x[(0:(k-1))-length(x)])


## UFD2: Chr 4-, 118-122000
chr = "4-"
sta = probeAnno$"4.-.start"
ind = probeAnno$"4.-.index"
uni = probeAnno$"4.-.unique"

sel  = (sta>=115000 & sta<=125000)
isel = ind[sel]
xsel = sta[sel]

qsel = (sta>=0) ##  & (sta<=300000)
qi   = ind[qsel]
qx   = sta[qsel]
xout = seq(min(qx), max(qx), by=8)

graphics.off()

if(!interact)
  sink("assessNorm.txt")

##
## colors
##
rgnf = range(refSig)
cols = colorRamp(c("red", "yellow", "blue"))((refSig[isel]-rgnf[1])/(rgnf[2]-rgnf[1])) / 256
## cols = colorRamp(brewer.pal(11, "Spectral"))((refSig[isel]-rgnf[1])/(rgnf[2]-rgnf[1])) / 256
cols = rgb(cols[,1], cols[,2], cols[,3])

ass = matrix(NA, nrow=5, ncol=length(rnaTypes))
colnames(ass)=rnaTypes

for(rt in rnaTypes) {
  fn = switch(rt,
    "seg-polyA-050525" = c("05_04_27_2xpolyA_NAP3.cel.gz",
      "05_04_26_2xpolyA_NAP2.cel.gz",
      "05_04_20_2xpolyA_NAP_2to1.cel.gz"),
    "seg-tot-050525"   = c("050409_totcDNA_14ug_no52.cel.gz",
      "030505_totcDNA_15ug_affy.cel.gz"),
    stop("Zapperlot"))

  ##
  ## plot
  ##
  if(makeFig){
    if(interact) {
      x11(width=8, height=12)
    } else {
      pdf(width=8, height=12, file=paste("assessNorm-", rt, ".pdf", sep=""))
    }
    
    py = calcY(isel, fn, rt)
    par(mfrow=c(ncol(py), 1))
    for(i in 1:ncol(py)) {
      sl  = !is.na(py[,i])
      plot(xsel[sl], py[sl,i], col=cols[sl], pch=20, xlab="coordinates", ylab="signal",
           main=colnames(py)[i])
    }
    if(!interact)
      dev.off()
  }
  
  ##
  ## quantitative assessment:
  ##
  py   = calcY(qi, fn, rt)
  max.lag = 100
  wh  = 50
  acr = matrix(NA, nrow=max.lag+1, ncol=ncol(py))
  colnames(acr) = colnames(py)
  
  for(i in 1:ncol(py)) {
    sl  = !is.na(py[,i])
    ap = approx(x=qx[sl], y=py[sl, i], xout=xout, ties=mean)
    if(FALSE){
      ac = acf(ap$y, lag.max=max.lag, na.action=na.pass, plot=FALSE)
      stopifnot(ac$lag[wh+1]==wh)
      acr[, i]   = ac$acf
    }
        
    dd = calcdd(ap$y)
    ## ee = calcdd(ap$y, k=1250)
    
    sig = diff(quantile(ap$y, probs=c(0.01, 0.99)))
    noi = median(dd)
    
    ass[i, rt] = sig / noi
    ## cat(rt, i, signif(c(sig, noi, ass[i, rt]), 3), "\n")
    stopifnot(!is.na(ass[i, rt]))
  }
  ## matplot(acr, type="l")
  
} ## for

cat("Signal-to-noise ratios:\n")
rownames(ass)=colnames(py)
print(signif(ass,3))
 
if(!interact)
  sink()
