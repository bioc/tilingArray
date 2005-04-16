library("tilingArray")
if(!exists("a"))load("a.rda")
if(!exists("x"))load("x.rda")
if(!exists("gff")) {
  load("probeAnno.rda")
  gff$Name = getAttributeField(gff$attributes, "Name")
}

source("/home/huber/madman/Rpacks/tilingArray/R/plotAlongChrom2.R")

fn  = "050209_mRNAx4_30min_re-hybe_RH6.cel.gz"

myPlot = function(y, vprow) {
  ## pushViewport(viewport(layout.pos.col=1, layout.pos.row=vprow))
  plotAlongChrom2(chr  = 14, coord = c(465000, 475000),
    y = y, probeAnno = probeAnno,
    gff  = gff, haveLegend=FALSE)
  ## popViewport()
  }

## pushViewport(viewport(layout=grid.layout(2, 1, height=c(1,1))))

## myPlot(log(exprs(a)[, fn], 2), 1)
myPlot(exprs(x)[, fn], 2)    

## popViewport()
