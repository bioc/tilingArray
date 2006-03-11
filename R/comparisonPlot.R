comparisonPlot = function(x, y, xscale=range(x), yscale) {
  
  myColorRamp = function(x, rg=range(x)) {
    cols = colorRamp(c("red", "yellow", "blue"))((x-rg[1])/(rg[2]-rg[1])) / 256
    return(rgb(cols[,1], cols[,2], cols[,3]))
  }
  cols = myColorRamp(x[[1]])

  n   = length(y)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(n+1, 1, height=c(rep(1,n), 0.3))))
  for(i in 1:n) { 
    pushViewport(viewport(layout=grid.layout(2, 2, height=c(.95, 0.05), width=c(0.1, 0.9)),
                            layout.pos.col=1, layout.pos.row=i, clip="off"))
    if(i<n) grid.lines(x=c(0,1), y=c(0,0), default.units="npc")
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=1, clip="off"))
    grid.text(letters[i], x=0, y=0.5, default.units="npc", hjust=0, vjust=0.5,
       gp=gpar(fontface="bold", fontsize=16))
    popViewport()
  
    pushViewport(dataViewport(xscale=xscale, yscale=yscale[,i],
                              layout.pos.col=2, layout.pos.row=1, clip="off"))
    grid.yaxis()
    if(i==n) grid.xaxis()
    popViewport()
  
    pushViewport(dataViewport(xscale=xscale, yscale=yscale[,i],
                              layout.pos.col=2, layout.pos.row=1, clip="on"))

    grid.points(px, y[[i]], pch=20, gp=gpar(col=cols))
    popViewport(2)
  }
  popViewport()

}
