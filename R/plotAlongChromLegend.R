##------------------------------------------------------------
## legend
##------------------------------------------------------------
plotAlongChromLegend = function(vpr, nr=2, 
    featureColorScheme=1, featureExclude=c("chromosome", "nucleotide_match", "insertion"),
    mainLegend, cexLegend=0.35, cexMain=1){
  
  endVP = FALSE
  # when this function is called on its own 
  # set up a viewport
  if(missing(vpr)) { 
     endVP=TRUE      
     vpr = newVP(main=mainLegend, dataPanelHeight=1, cexMain=cexMain) # newVP sets up a new viewport
  }
  formatRow = function(featColsOneRow, row) {
    ## print(featColsOneRow)
    strWid   = convertWidth(stringWidth(rownames(featColsOneRow)), "npc", valueOnly=TRUE)
    n        = length(strWid)
    inbetWid = 0.2*min(strWid)
    totWid   = sum(strWid)+(n-1)*inbetWid
    x        = c(0, cumsum(strWid[-n])) + (0:(n-1))*inbetWid 
    y        = numeric(length(x))

    x      = x/totWid
    strWid = strWid/totWid
    grid.rect(x = x, width = strWid, 
              y = unit(row, "native"), height = unit(1, "native")- unit(1, "mm"), 
              just  = c("left", "center"), default.units="npc",
              gp    = do.call("gpar", featColsOneRow))
    
    grid.text(label = rownames(featColsOneRow),
              x = unit(x + strWid/2, "native"), y = unit(row, "native"),
              just  = c("center", "center"), gp=gpar(cex=cexLegend))
  } 

  featCols = featureColors(featureColorScheme)
  featCols = featCols[ !(rownames(featCols) %in% featureExclude), ]

  pushViewport(viewport(layout.pos.col=1, layout.pos.row=vpr, yscale=c(0.5, nr+0.5)))

  i = 1:nrow(featCols)
  for(r in 1:nr)
    formatRow(featCols[ceiling(i/nrow(featCols)*nr-1e-10)==r, ], row=nr-r+1)
  
  popViewport()

  if(endVP)
     popViewport(2)
}

# this function sets up a new viewport.  It is used by plotAlongChromLegend, 
# plotSegmentationHeatmap and plotSegmentationDots when they are called as 
# stand-alone functions (ie when vpr is not specified)

newVP <- function(main, cexMain=1, dataPanelHeight=1, vpHeight=0.7, titleOffSet=0) {
     if(!missing(main)) {
        vpr = c("title"=0.1, "data"=dataPanelHeight)
        pushViewport(viewport(width=0.85, height=vpHeight)) ## plot margin
        pushViewport(viewport(layout=grid.layout(length(vpr), 1, height=vpr)))  
        pushViewport(viewport(layout.pos.col=1, layout.pos.row=which(names(vpr)=="title")))
        grid.text(label=main, x=0.5, y=1.1+titleOffSet, just="centre", gp=gpar(cex=cexMain))  
        popViewport()
        vpr = which(names(vpr)=="data")
     } else {
        vpr = c("data"=dataPanelHeight)
        pushViewport(viewport(width=0.85, height=vpHeight)) ## plot margin
        pushViewport(viewport(layout=grid.layout(length(vpr), 1, height=vpr)))
     }
  vpr
  }
