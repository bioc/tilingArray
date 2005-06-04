options(error=recover, warn=0)
library("tilingArray")
library("geneplotter")

## source("colorRamp.R")  ## can go with R 2.1
source("~/madman/Rpacks/tilingArray/R/plotAlongChrom2.R")

source("scripts/readSegments.R")
doSave = TRUE

so = get("seg-polyA-050525")


density = 100
width   = 10
height  = 12
graphics.off()
wkb     = 7 ## width in kilobases
ylim    = c(-4.5, 3.5)

if(doSave) {
  fn = "Figures/fig1.pdf"
  pdf(file=fn, width=width, height=height)
} else {
  x11(width=width, height=height)
}

grid.newpage()

## push 1: first level of layout
pushViewport(viewport(layout=grid.layout(3, 2, height=c(1, 0.1, 2), width=c(0.02, 1))))

##
## 13:541-555: the segmentation nicely captures two spliced transcripts
##
## push 2: top panel
pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
plotAlongChrom2(main="a)", chr=13, coord = (541+c(0, 2*wkb))*1e3, ylim=ylim, segObj=so, 
                probeAnno = probeAnno, gff=gff, haveLegend=FALSE)
## pop 2
popViewport()

## push 2: second level of layout
pushViewport(viewport(layout.pos.col=2, layout.pos.row=3,
      layout=grid.layout(3, 3, height=c(1, 0.1, 1), width=c(1, 0.1, 1))))

##
## 1:41-47 ACS1: novel architecture
##
## push 3: middle left
pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
plotAlongChrom2(main="b)", chr=1, coord = (41+c(0, wkb))*1e3, ylim=ylim, segObj=so, 
                probeAnno = probeAnno, gff=gff, haveLegend=FALSE)
## pop3
popViewport()

##
## 11:62.5-69.5: MNN4, novel architecture
##
## push 3: middle right
pushViewport(viewport(layout.pos.col=3, layout.pos.row=1))
plotAlongChrom2(main="c)", chr=11, coord = (62.5+c(0, wkb))*1e3, ylim=ylim, segObj=so, 
                probeAnno = probeAnno, gff=gff, haveLegend=FALSE)
## pop3
popViewport()

##
## 2:360-367: novel isolated
##
## push 3: middle right
pushViewport(viewport(layout.pos.col=1, layout.pos.row=3))
plotAlongChrom2(main="d)", chr=2, coord = (360+c(0, wkb))*1e3, ylim=ylim, segObj=so, 
                probeAnno = probeAnno, gff=gff, haveLegend=FALSE)
## pop3
popViewport()

##
## 6:87-94: novel antisense
##
## push 3: middle right
pushViewport(viewport(layout.pos.col=3, layout.pos.row=3))
plotAlongChrom2(main="e)", chr=6, coord = (87+c(0, wkb))*1e3, ylim=ylim, segObj=so, 
                probeAnno = probeAnno, gff=gff, haveLegend=FALSE)
## plotAlongChrom2(chr=1, coord = c(141.2+144.5)*1e3, ylim=ylim, segObj=so, 
##                probeAnno = probeAnno, gff=gff, haveLegend=FALSE)
## pop3
popViewport()

for(j in c(1,3)) {
  pushViewport(viewport(layout.pos.row=j, layout.pos.col=2))
  grid.lines(x=0.5, y=c(0,1), gp=gpar(col="grey"))
  popViewport()
  pushViewport(viewport(layout.pos.row=2, layout.pos.col=j))
  grid.lines(x=c(0,1), y=0.5, gp=gpar(col="grey"))
  popViewport()
}

## pop 2,1
popViewport(2)

if(doSave) {
  dev.off()
  cmd = paste("convert -density 180", fn, "-compress RLE", sub(".pdf", ".tiff", fn))
  system(cmd)
}
