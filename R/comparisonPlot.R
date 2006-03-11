comparisonPlot = function(x, y, xscale=range(x), yscale, anno, ticks) {
  
  myColorRamp = function(d, rg=range(d)) {
    cols = colorRamp(c("red", "yellow", "blue"))((d-rg[1])/(rg[2]-rg[1])) / 256
    return(rgb(cols[,1], cols[,2], cols[,3]))
  }
  cols = myColorRamp(y[[1]])
  n = length(y)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(n+2, 1, height=c(rep(1,n), 0.2, 0.2))))
  for(i in 1:n) { 
    pushViewport(viewport(layout=grid.layout(2, 2, height=c(.95, 0.05), width=c(0.1, 0.9)),
                            layout.pos.col=1, layout.pos.row=i, clip="off"))
    ## if(i<n) grid.lines(x=c(0,1), y=c(0,0), default.units="npc")
    if(i%%2==1) grid.rect(x=0.5, y=0.5, width=1, height=1, default.units="npc",
            gp=gpar(col="#f0f0f0", fill="#f0f0f0"))
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=1, clip="off"))
    grid.text(letters[i], x=0, y=0.5, default.units="npc", hjust=0, vjust=0.5,
       gp=gpar(fontface="bold", fontsize=16))
    popViewport()
  
    pushViewport(dataViewport(xscale=xscale, yscale=yscale[,i],
                              layout.pos.col=2, layout.pos.row=1, clip="off"))
    grid.yaxis()
    popViewport()
  
    pushViewport(dataViewport(xscale=xscale, yscale=yscale[,i],
                              layout.pos.col=2, layout.pos.row=1, clip="on"))

    grid.points(x, y[[i]], pch=".", gp=gpar(col=cols))
    popViewport(2)
  }
  pushViewport(viewport(layout=grid.layout(2, 2, height=c(.95, 0.05), width=c(0.1, 0.9)),
                            layout.pos.col=1, layout.pos.row=n+1, clip="off"))
  pushViewport(dataViewport(xscale=xscale, yscale=c(0,1),
                              layout.pos.col=2, layout.pos.row=1, clip="off"))
  grid.lines(xscale, y=c(0.5,0.5), "native")
  for (i in 1:nrow(anno)) {
    grid.rect(x=anno$start[i], y=0.5, width=anno$end[i]-anno$start[i]+1, height=0.7, default.units="native",
              just  = c("left", "center"), gp=gpar(col="black", fill="#d0d0d0"))
    grid.text(x=(anno$start[i]+anno$end[i])/2, y=0.5, label=anno$name[i], default.units = "native")
  }
  grid.segments(x0=ticks*1e3, x1=ticks*1e3, y0=0.4, y1=0.6, default.units="native")
  grid.text(x=ticks*1e3, y=0.3, label=sprintf("%dkB", ticks), just=c("center", "top"), default.units="native")
  popViewport(3)
}

