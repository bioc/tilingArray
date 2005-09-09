library("tilingArray")
source("setScriptsDir.R")

graphics.off()
options(error=recover, warn=2)

interact = TRUE
rnaTypes = c("seg-polyA-050811", "seg-tot-050811")
outfile  = "lowTranscribedRegions"
doNotLoadSegScore = TRUE

source(scriptsDir("readSegments.R")) 

if(!interact){
  sink(paste(outfile, ".txt", sep=""))
  cat("Made on", date(), "\n\n")
}

movingAverage = function(x, y, win, by=4, nmin=4) {
  myMean = function(x) {
    res=as.numeric(NA)
    if(length(x)>=nmin)
      res=sum(x)/length(x)
    return(res)
  }
        
  xout = seq(min(x)+win/2, max(x), by=win/by)
  yout = numeric(length(xout))
  mx   = max(xout)
  
  for(j in 1:by) {
    isel = seq(j, length(xout), by=by)
    ct   = cut(x, breaks=c(xout[isel]-win/2, mx), include.lowest=TRUE)
    yout[isel] = tapply(y, ct, myMean)
    ##browser()
  }
  
  return(list(xout=xout,yout=yout))
}


####
windowWidth = 100
chrs = 1:16

for(rt in rnaTypes) {
  mavs = list()
  for(chr in chrs) {
    for(strand in c("+", "-")) {
       dat = get(paste(chr, strand, "dat", sep="."), envir=get(rt))

       x  = (dat$start+dat$end)/2
       y  = rowMeans(dat$y)

       ma = movingAverage(x, y, win=windowWidth)
       mavs[[paste(chr, strand, sep=".")]] = ma

       plot(x, y, pch=".", main=paste(chr, strand, sep="."))
       points(ma$xout, ma$yout, pch=".", col="orange")

    }
  }
}

if(!interact)
  sink()
