grid.image = function(x, y, xlim, uniq=uniq,
  colRamp=colorRamp(c("blue", "yellow")),
  transformation=function(z) { res=rank(z)/length(z); dim(res)=dim(z); return(z) },
  width=6) {
  
  pushViewport(dataViewport(xscale=xlim, yscale=c(.5, ncol(y)+.5)))  

  mcol = (colRamp(transformation(z)))/256
  mcol[uniq!=0,] = rep(255, 3)/256 # colour of non-unique probes - white was grey: 190
  
  col = rgb(mcol[,1], mcol[,2], mcol[,3])
  grid.rect(rep(x, ncol(y)), rep(1:ncol(y), each=nrow(y)), 
             width=width, height=1, default.units="native",
             gp=gpar(col=col, fill=col) )
  popViewport()
}
