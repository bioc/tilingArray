##------------------------------------------------------------- 
##   plotSegmentationHeatmap
##------------------------------------------------------------- 
plotSegmentationHeatmap = function(dat, xlim, ylab, rowNames, 
  chr=1, strand="+", vpr, colors, colHeatmap=colorRamp(brewer.pal(9, "YlGnBu")),
  main,...) {

  endVP = FALSE
  if(missing(vpr)) {
     endVP=TRUE
     vpr = newVP(main=main, dataPanelHeight=1, vpHeight=0.95, titleOffSet=-0.9)
  }

  if(min(dat$y, na.rm=TRUE)<0 | max(dat$y, na.rm=TRUE)>1) {
    ythresh = quantile(dat$y, c(0.05,0.95), na.rm=TRUE)
    ymod = dat$y[,,drop=FALSE]
    ymod[ymod<ythresh[1]] = ythresh[1]
    ymod[ymod>ythresh[2]] = ythresh[2]
    dat$y = (ymod - min(ymod, na.rm=TRUE))/(max(ymod, na.rm=TRUE)-min(ymod, na.rm=TRUE))
  }

  if(missing(xlim)) {
    xlim=range(dat$x, na.rm=TRUE)
  } else {
    sel = (dat$x>=xlim[1])&(dat$x<=xlim[2])
    dat$x = dat$x[sel]
    dat$y = dat$y[sel,, drop=FALSE ]
    dat$flag = dat$flag[sel]
  }

  defaultColors = c("+" = "#00441b", "-" = "#081d58", "duplicated" = "grey",
                    "cp" = "#555555", "ci" = "#777777", "highlight" = "red", "threshold" = "grey")

  if(!missing(colors)) {
    mt = match(names(colors), names(defaultColors))
    if(any(is.na(mt)))
      stop(paste("Cannot use color specification for", names(colors)[is.na(mt)]))
    defaultColors[mt] = colors 
  }
  colors = defaultColors

  ord = order(dat$x)
  dat$x = dat$x[ord]   ## sort by x-coordinates to simplify smoothing
  dat$y = dat$y[ord,, drop=FALSE]
  dat$flag = dat$flag[ord]
  
  ## Use two viewports for different clipping behavior
  ylim = c(-1, 2+ncol(dat$y))
  pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="off",
    layout.pos.col=1, layout.pos.row=vpr))

  if(missing(rowNames))
    rowNames = colnames(dat$y)
  if(!is.null(rowNames))
    grid.yaxis( (1:ncol(dat$y)), rowNames, gp=gpar(cex=0.5))

  if(!missing(ylab))
    grid.text(ylab, x=-0.075, y=0.5, rot=90)

  pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="on",
    layout.pos.col=1, layout.pos.row=vpr))

  ord  = c(which(dat$flag!=0), which(dat$flag==0))
  colo = ifelse(dat$flag[ord]==0, colors[strand], colors["duplicated"])

  grid.image(dat$x, 1:ncol(dat$y), z=dat$y, xlim=xlim, uniq=dat$flag, colRamp=colHeatmap)

  ## segment boundaries
  if(!is.null(dat$estimate)) {
    grid.segments(x0 = unit(dat$estimate, "native"),
                  x1 = unit(dat$estimate, "native"),
                  y0 = unit(0.1, "npc"),
                  y1 = unit(0.9, "npc"),
                  gp = gpar(col=colors["cp"]))
  }

  popViewport(2)  

  if(endVP)
     popViewport(2)

} ## end of plotSegmentationHeatmap
