options(error=recover, warn=0)

library("tilingArray")

## source("colorRamp.R")  ## can go with R 2.1
source("~/madman/Rpacks/tilingArray/R/plotAlongChrom2.R")


#### Generic plot
graphics.off();
## X11(width=15, height=8); grid.newpage()
pdf(file="test-utr.pdf", width=15, height=8)

if(!TRUE) {
  ## with segRes environment
  source("scripts/readSegments.R")
  rt = "polyA"
  plotAlongChrom2(which(gff$seqname[w]==chrSeqname), coord = c(gff$start[w]-1e4, gff$end[w]+1e4),
                nrBasesPerSeg=1500, segRes = get(rt),
                ## segScore = get("segScore", e), 
                gff = gff, highlight= list(coord=c(142621, 143365),strand="+"))
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
      plotAlongChrom2(2, coord = c(start[i], end[i])*1e3,
                  y = zz[,j], probeAnno = probeAnno,
                  isDirectHybe = FALSE, 
                  gff = gff)
      popViewport()
    }
  popViewport()
  
}

dev.off()
