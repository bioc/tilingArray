doSave = TRUE
interact = FALSE

options(error=recover, warn=0)
library("tilingArray")
library("geneplotter")

source("~/madman/Rpacks/tilingArray/R/plotAlongChrom2.R")

rnaTypes  = c("seg-polyA-050525", "seg-tot-050525")
source("scripts/readSegments.R")
source("scripts/calcThreshold.R") 



so = get("seg-polyA-050525")

density = 100
width   = 10
height  = 12
graphics.off()
wkb     = 6 ## width in kilobases
ylim    = c(-4.5, 3.5)

if(doSave) {
  fn = "Figures/fig1.pdf"
  pdf(file=fn, width=width, height=height)
} else {
  x11(width=width, height=height)
}

grid.newpage()

## push 1: first level of layout
pushViewport(viewport(layout=grid.layout(2, 1, height=c(1.5, 2), width=1)))

##
## 14:380-480: overview
##
## push 2: top panel
pushViewport(viewport(layout.pos.col=1, layout.pos.row=1,
     layout=grid.layout(1, 2, height=1, width=c(0.01, 1))))
pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))

plotAlongChrom2(main="a)", chr=14, coord= c(410,490)*1e3,
                ylim=ylim, segObj=so, colors=c(cp="#d0d0d0"), 
                probeAnno = probeAnno, gff=gff,
                haveNames=FALSE, haveLegend=FALSE, pointSize=unit(0.1, "mm"))
## pop 2
popViewport(2)

## -------------------------------------------------------------
## push 2: second level of layout
pushViewport(viewport(layout.pos.col=1, layout.pos.row=2,
      layout=grid.layout(2, 6, height=c(1, 1), width=c(0.1, 1, 0.07, 1, 0.07, 1))))

##
## 13:548.6-554.6 splicing RPS16A, RPL13B
##
## push 3: middle left
pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
plotAlongChrom2(main="b)", chr=13, coord = (548.6+c(0, wkb))*1e3, ylim=ylim, segObj=so, 
                probeAnno = probeAnno, gff=gff, haveLegend=FALSE)
## pop3
popViewport()

##
## 11:63-69: MNN4, novel architecture
##
## push 3: middle right
pushViewport(viewport(layout.pos.col=4, layout.pos.row=1))
plotAlongChrom2(main="c)", chr=11, coord = (63+c(0, wkb))*1e3, ylim=ylim, segObj=so, 
                probeAnno = probeAnno, gff=gff, haveLegend=FALSE)
## pop3
popViewport()

##
## 2:360.5-366.5: novel isolated
##
## push 3: bottom left
pushViewport(viewport(layout.pos.col=6, layout.pos.row=1))
plotAlongChrom2(main="d)", chr=2, coord = (360.5+c(0, wkb))*1e3, ylim=ylim, segObj=so, 
                probeAnno = probeAnno, gff=gff, haveLegend=FALSE)
## pop3
popViewport()

##
## 9:221-227: novel antisense SPO22
##
## push 3: bottom right
pushViewport(viewport(layout.pos.col=2, layout.pos.row=2))
plotAlongChrom2(main="e)", chr=9, coord = (221+c(0, wkb))*1e3, ylim=ylim, segObj=so, 
                probeAnno = probeAnno, gff=gff, haveLegend=FALSE)
## plotAlongChrom2(chr=1, coord = c(141.2+144.5)*1e3, ylim=ylim, segObj=so, 
##                probeAnno = probeAnno, gff=gff, haveLegend=FALSE)
## pop3
popViewport()

#for(j in c(1,3)) {
#  pushViewport(viewport(layout.pos.row=j, layout.pos.col=2))
#  grid.lines(x=0.5, y=c(0,1), gp=gpar(col="grey"))
#  popViewport()
#  pushViewport(viewport(layout.pos.row=2, layout.pos.col=j))
#  grid.lines(x=c(0,1), y=0.5, gp=gpar(col="grey"))
#  popViewport()
#}

## pop 2,1
popViewport(2)

if(doSave) {
  dev.off()
  cmd = paste("convert -density 180", fn, "-compress RLE", sub(".pdf", ".tiff", fn))
  cat(cmd, "\n")
  system(cmd)
}
