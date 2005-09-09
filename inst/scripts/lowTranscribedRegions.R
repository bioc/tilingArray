
nrSelect  = 250

library("tilingArray")
library("arrayMagic")

source("setScriptsDir.R")

graphics.off()
options(error=recover, warn=2)

interact = TRUE
rnaTypes = c("seg-polyA-050811", "seg-tot-050811")[1]


stopifnot(length(rnaTypes)==1)
rt = rnaTypes
outfile  = file.path(rt, "lowTranscribedRegions")
doNotLoadSegScore = TRUE

source(scriptsDir("readSegments.R")) 

movingAverage = function(x1, y1, x2, y2, windowWidth, windowStep, nmin=4) {
  stopifnot(all(diff(x1)>=0), all(diff(x2)>=0))
  rx   = range(c(x1, x2))
  xout = seq(rx[1]+windowWidth, rx[2]-windowWidth, windowStep)
  xs   = xout-windowWidth/2
  xe   = xout+windowWidth/2
  
  yout = rep(as.numeric(NA), length(xout))
  i1 = i2 = j1 = j2 = 1
  
  for(k in seq(along=xout)) {
    while(x1[i1] < xs[k])
      i1=i1+1 
    while(x2[i2] < xs[k])
      i2=i2+1
    while(x1[j1] <= xe[k])
      j1=j1+1 
    while(x2[j2] <= xe[k])
      j2=j2+1
    ##cat(k, "\n",
    ##    i1, j1, xs[k], "<=", x1[i1], "<=", x1[j1-1], "<=", xe[k], "\n",
    ##    i2, j2, xs[k], "<=", x2[i2], "<=", x2[j2-1], "<=", xe[k], "\n")
    if( (j1-i1 >= nmin) && (j2-i2 >= nmin) )
      yout[k] = mean(c(y1[i1:(j1-1)], y2[i2:(j2-1)]))
    ## browser()
  }
  return(list(xout=xout,yout=yout))
}


####
windowWidth = 100
windowStep  = 25
chrs = 1:16

if(!exists("mavs")) {
mavs = vector(mode="list", length=length(chrs))
for(chr in chrs) {
  dat1 = get(paste(chr, "+.dat", sep="."), envir=get(rt))
  dat2 = get(paste(chr, "-.dat", sep="."), envir=get(rt))
  
  x1  = (dat1$start+dat1$end)/2
  x2  = (dat2$start+dat2$end)/2
  y1  = rowMeans(dat1$y)
  y2  = rowMeans(dat2$y)
  
  mavs[[chr]] = movingAverage(x1, y1, x2, y2, windowWidth, windowStep)
  
  plot(x1, y1, pch=".", col="grey", main=paste(chr, sep="."))
  points(x2, y2, pch=".", col="grey")
  points(ma$xout, ma$yout, pch=".", col="orange")
}
}

allValues = unlist(lapply(mavs, "[[", "yout"))
thresh    = sort(allValues)[nrSelect]

res = data.frame(chr=integer(nrSelect),
  midpoint=integer(nrSelect),
  plot=I(character(nrSelect)))


alongChromWidth = 15e3
alongChromStep  =  5e3

## this function maps a chromosome number and start and end coordinates
## to a file name
mapCoord2Plot = function (chr, start, end) {
  mid  = (start+end-alongChromWidth)/2
  mid[mid<0]=0
  pst  = as.integer(alongChromStep/1e3*ceiling(mid/alongChromStep))
  sprintf("%02d_%04d", as.integer(chr), pst)
}

j = 0
for(chr in chrs) {
  wh = which(mavs[[chr]]$yout <= thresh)
  if(length(wh)>0){
    mids = mavs[[chr]]$xout[wh]
    idx = j + (1:length(mids))
    
    res[idx, "chr"] = chr
    res[idx, "midpoint"] = mids
    hyper = mapCoord2Plot(chr, mids, mids)
    res[idx, "plot"] = sprintf("<a href=\"viz/%s.jpg\" target=\"new\">%s</a>",
         hyper, hyper)
    j = j + length(mids)
  }
}

write.htmltable(res, outfile, title=rt, formatNumeric=function(x) paste(as.integer(x)))
