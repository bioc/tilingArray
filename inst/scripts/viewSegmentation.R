##
## This script presents an example for calling the plotAlongChrom function
##

options(error=recover, warn=0)
library("tilingArray")

source("setScriptsDir.R")


source(functionsDir("plotAlongChrom.R"))
source(scriptsDir("readSegments.R"))

graphics.off()

out = c("x11", "pdf")[1]
switch(out,
       x11 = {
         X11(width=15, height=8)
         grid.newpage()
       },
       pdf = {
         pdf(file="viewSeg.pdf", width=15, height=8)
       },
       stop("Sapperlot"))



if(TRUE) {
  ## with segRes environment
  rnaTypes = rt = "seg-polyA-050909"
  source("scriptsd/readSegments.R")
  plotAlongChrom(chr=1, coord = 1000*c(30, 130),
                segObj = get(rt),
                gff = gff, isDirect=FALSE)
} else {
  ## if(!exists("a"))load("a.rda")
  if(!exists("probeAnno"))load("probeAnno.rda")

  if(!exists("xn1")) {
    load("seg-polyA-050525/xn.rda")
    xn1=xn
  }
  if(!exists("xn2")) {
    load("seg-tot-050525/xn.rda")
    xn2=xn
  }

  zz = cbind(rowMeans(exprs(xn1)), rowMeans(exprs(xn2)))
  start = c(113.5, 115.3)
  end   = c(113.8, 115.6)

  pushViewport(viewport(layout=grid.layout(2, 2)))
  for(i in 1:2)
    for(j in 1:2) {
      pushViewport(viewport(layout.pos.col=i, layout.pos.row=j))
      plotAlongChrom(2, coord = c(start[i], end[i])*1e3,
                  y = zz[,j], probeAnno = probeAnno,
                  isDirectHybe = FALSE, 
                  gff = gff)
      popViewport()
    }
  popViewport()
  
}

if(out %in% "pdf")
  dev.off()
