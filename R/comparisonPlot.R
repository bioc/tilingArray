comparisonPlot = function(x, y, xscale=range(x), yscale, anno, ticks, pch=20, cex=1, bgcol="#f2f2f2") {
  
  myColors = function(z){
    brewer.pal(6, "Set1")[c(1,5,3,2)][as.integer(cut(z, quantile(z, c(0, 0.05, 0.35, 0.65, 1)), include.lowest=TRUE))]
  }  
  cols = myColors(rank(y[[1]]))
  
  n = length(y)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(n+2, 1, heights=c(rep(1,n), 0.2, 0.2))))
  for(i in 1:n) { 
    pushViewport(viewport(layout=grid.layout(2, 2, heights=c(.95, 0.05), widths=c(0.1, 0.9)),
                            layout.pos.col=1, layout.pos.row=i, clip="off"))
    if(i%%2==1) grid.rect(x=0.5, y=0.5, width=1, height=1, default.units="npc",
            gp=gpar(col=bgcol, fill=bgcol))
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

    grid.points(x, y[[i]], pch=pch, gp=gpar(col=cols))
    popViewport(2)
  }
  pushViewport(viewport(layout=grid.layout(2, 2, heights=c(.95, 0.05), widths=c(0.1, 0.9)),
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

