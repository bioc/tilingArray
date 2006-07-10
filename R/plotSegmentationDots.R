##-------------------------------------------------------------
##   plot Segmentation with Dots
##-------------------------------------------------------------
plotSegmentationDots = function(dat, xlim, ylim, ylab, threshold=NA, 
  chr=1, strand="+", vpr, colors, main, pointSize=unit(0.6, "mm"),
  showConfidenceIntervals=TRUE, sepPlots=FALSE, cexAxisLabel=1, cexAxis=1,...) {

  endVP = FALSE
  if(missing(vpr)) {
     endVP=TRUE
     vpr = newVP(main=main, dataPanelHeight=1, vpHeight=0.95, titleOffSet=-0.9)
  }

  if(is.matrix(dat$y) & !sepPlots)
    dat$y = rowMeans(dat$y) ##  if >1 samples, take mean over samples
  stopifnot(length(dat$y)==length(dat$x), length(dat$flag)==length(dat$x))
  
  xorg  = dat$x
  if(missing(xlim)) {
    xlim=range(dat$x, na.rm=TRUE)
  } else {
    sel  = (dat$x>=xlim[1])&(dat$x<=xlim[2])
    dat$x = dat$x[sel]
    dat$y = dat$y[sel]
    dat$flag = dat$flag[sel]
  }
  
  if(!is.na(threshold)) {
    dat$y = dat$y-threshold
    if(!missing(ylim))
      ylim = ylim-threshold
  }
 
  if(missing(ylim))
    ylim = quantile(dat$y, c(0,1), na.rm=TRUE)

  ## the expression data. use two viewports for different clipping behavior
  if(!missing(ylab)) {
    pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="off",
      layout.pos.col=1, layout.pos.row=vpr))
    grid.yaxis(gp=gpar(cex=cexAxis),...)
    grid.text(ylab, x=-0.075, y=0.5, rot=90, gp=gpar(cex=cexAxisLabel),...)
    pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="on",
      layout.pos.col=1, layout.pos.row=vpr))
    }
  else {
    pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="off",
      layout.pos.col=1, layout.pos.row=vpr))
    grid.yaxis(gp=gpar(cex=cexAxis))
    pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="on",
      layout.pos.col=1, layout.pos.row=vpr))
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

  ord  = c(which(dat$flag!=0), which(dat$flag==0))
  colo = ifelse(dat$flag[ord]==0, colors[strand], colors["duplicated"])

  if(!is.na(threshold))
    grid.lines(y=unit(0, "native"), gp=gpar(col=colors["threshold"]))

  ## segment boundaries
  sel = ((xorg[dat$estimate]>=xlim[1]) & (xorg[dat$estimate]<=xlim[2]))
  mySeg = function(j, what)
    grid.segments(x0 = unit(j, "native"), x1 = unit(j, "native"),
                  y0 = unit(0.1, "npc"),  y1 = unit(0.9, "npc"),
                  gp = gpar(col=colors[what], lty=c(cp=1, ci=2)[what]))
    
  if(!is.null(dat$estimate) & sum(sel)>0) mySeg(xorg[dat$estimate][sel], "cp")
  if(showConfidenceIntervals & sum(sel)>0) {
    if(!is.null(dat$upper))    mySeg(xorg[dat$upper][sel], "ci")
    if(!is.null(dat$lower))    mySeg(xorg[dat$lower][sel], "ci")
  }
  
  grid.points(dat$x[ord], dat$y[ord], pch=20, size=pointSize, gp=gpar(col=colo))
  popViewport(2)

  if(endVP)
     popViewport(2)

} ## plotSegmentationDots
