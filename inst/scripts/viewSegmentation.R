##
## This script presents an example for calling the plotAlongChrom function
##
## It is also used for the plotAlongChrom figure in the CompStatViz Paper
##
options(error=recover, warn=0)
library("tilingArray")
source("setScriptsDir.R")

## source(functionsDir("plotAlongChrom.R"))
## source(functionsDir("grid.image.R"))
options(error=recover)

graphics.off()

what = c("dotsSeg", "dotsUnseg", "heatmap")[3]
name = sprintf("fig_tiling_%s", what)

for(dev in c("png", "pdf")) {
  switch(dev,
         X11 = {
           X11(width=15, height=8)
           grid.newpage()
         },
         pdf = pdf(file=sprintf("%s.pdf", name), width=8, height=4.7),
         png = png(file=sprintf("%s.png", name), width=1024, height=768),
         stop("Sapperlot"))
  
  
  ## dots plot with segmentation
  if(what=="dotsSeg") {
    rnaTypes = rt = "seg-polyA-050909"
    source(scriptsDir("readSegments.R"))
    source(scriptsDir("calcThreshold.R"))
    plotAlongChrom(chr=1, coord = 1000*c(76.2, 92.8),
                   segObj = get(rt),
                   gff = gff)
  }
  
  ## dots plot without segmentation, just normalized intensities
  if(what=="dotsUnseg"){
    ## if(!exists("a"))load("a.rda")
    if(!exists("probeAnno"))load("probeAnno.rda")
    
    if(!exists("xn1")) {
      load("seg-polyA-050909/xn.rda")
      xn1=xn
    }
    if(!exists("xn2")) {
      load("seg-tot-050909/xn.rda")
      xn2=xn
    }
    
    zz = cbind(rowMeans(exprs(xn1)), rowMeans(exprs(xn2)))
    start = c(113.5, 115.3)
    end   = c(113.8, 115.6)
    
    pushViewport(viewport(layout=grid.layout(2, 2)))
    for(i in 1:2)
      for(j in 1:2) {
        pushViewport(viewport(layout.pos.col=i, layout.pos.row=j))
        plotAlongChrom(chr=2, coord = c(start[i], end[i])*1e3,
                       y = zz[,j,drop=FALSE], probeAnno = probeAnno,
                       gff = gff)
        popViewport()
      }
    popViewport()
    
  }
  
  ## heatmap plot
  if(what=="heatmap"){
    if(!exists("ex")) {
      load("/ebi/research/huber/Projects/tilingCycle/xn.rda")
      ex = exprs(xn)
      colnames(ex) = xn$SampleID
      time = as.numeric(colnames(ex))
      stopifnot(!any(is.na(time)))
      ex = ex[, order(time)]
    }
    if(!exists("probeAnno"))
      load("probeAnno.rda")
    
    smoothRank = function(y) {
      stopifnot(nrow(y)>=5)
      idx = 3:(nrow(y)-2)
      res = matrix(as.numeric(NA), nrow=nrow(y), ncol=ncol(y))
      res[idx, ] = y[idx,] + y[idx+1,] + y[idx+2,] + y[idx-1,] + y[idx-2,]
      return(rank(res) / length(res))
    }
    
    justRank = function(y) {
      return(rank(y) / length(y))
    }
    
    plotAlongChrom(y = ex, probeAnno = probeAnno, gff=gff,
                   what = "heatmap",
                   transformation = smoothRank,
                   chr= 1, coord = c(60, 66.5)*1e3)
    
  }
  if(dev %in% c("pdf", "png"))
    dev.off()
} ## for    







          
