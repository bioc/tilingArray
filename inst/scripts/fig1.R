doSave = TRUE
interact = !doSave

options(error=recover, warn=0)
library("tilingArray")
library("geneplotter")

source("~/madman/Rpacks/tilingArray/R/plotAlongChrom.R")

rnaTypes  = c("seg-polyA-050909")
source("scripts/readSegments.R")
source("scripts/calcThreshold.R") 

stopifnot(length(rnaTypes)==1)
so = get(rnaTypes)

density = 100
width   = 11
height  = 6
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
 

if(!exists("myGff"))
  myGff = gff[ gff$Name!="tR(UCU)E", ]

## -------------------------------------------------------------
## layout
dx = 0.18
dy = 0.02
pushViewport(viewport(x=0.01, width=0.97, height=0.97, just=c("left", "center"),
                      layout=grid.layout(3, 8,
                        height=c(1, dy, 1),
                        width =c(dx, 1, dx, 1, dx, 1, dx, 1))))

myPlot = function(row, col, chr, coord) {
  pushViewport(viewport(layout.pos.row=row, layout.pos.col=col))
  grid.rect(x=-0.1, width=1.15, y=0, height=1, just=c("left", "bottom"), default.units="npc")
  plotAlongChrom(chr=chr, coord = coord, ylim=ylim, segObj=so, 
                 probeAnno = probeAnno, gff=myGff,
                 noTypeLabel = c("CDS", "uORF", "binding_site", "TF_binding_site"),
                 haveLegend=FALSE)  
  popViewport()
}


##
## A) 13:550k splicing RPS16A, RPL13B
myPlot(1, 2, chr=13, coord = c(550044, 553360))

## B) GCN4
myPlot(1, 4, chr=5, coord = c(138660, 141880))

## C) MET7, novel architecture
myPlot(1, 6, chr=15, coord = c(784700, 790000))

## D) overlapping transcripts
myPlot(1, 8, chr=14, coord = c(342200, 347545))

## E) SER3
myPlot(3, 2, chr=5, coord = c(321900, 326100))

## F) 2:360.5-366.5: novel isolated
myPlot(3, 4, chr=2, coord = c(360500, 365970))

## G) 9:221-227: novel antisense SPO22
myPlot(3, 6, chr=9, coord = c(221000, 226500))

## separator line
## for(j in 1:8) {
##  pushViewport(viewport(layout.pos.col=j, layout.pos.row=2))
##  grid.lines(c(0,1), c(1,1), default.units = "npc", gp=gpar(col="#a0a0a0", lty=1, lwd=1))
##  popViewport()
##}

## legend
fc = featureColors(1)[c("CDS", "CDS_dubious", "uORF", "ncRNA", "TF_binding_site"), ]
pushViewport(dataViewport(xscale=c(0,1), yscale=c(-7,nrow(fc)+1), layout.pos.col=8, layout.pos.row=3))
h  = nrow(fc):1 
grid.rect(x = 0, width = 0.2, y = h, height = unit(1, "native")- unit(2, "mm"), 
            just  = c("left", "center"), default.units="native",
            gp    = do.call("gpar", fc))
grid.text(label = rownames(fc), x = 0.22, y = h,
            just  = c("left", "center"), default.units="native",
            gp=gpar(cex=.7))
popViewport()


popViewport()

if(doSave) {
  dev.off()
  cmd = paste("convert -density 180", fn, "-compress RLE", sub(".pdf", ".tiff", fn))
  cat(cmd, "\n")
  if(TRUE) {
    cat("PLEASE RUN!\n")
  } else {
    system(cmd)
  }
}
