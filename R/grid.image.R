grid.image = function(x, y, z, xlim, uniq=uniq,
  colRamp=colorRamp(brewer.pal(9, "YlGnBu")),
  width=22) {
  
  mcol = (colRamp(z))/256
  mcol[uniq!=0,] = 1
  indna <- is.na(mcol[,1])
  mcol[indna,] = 1

  col = rgb(mcol[,1], mcol[,2], mcol[,3])
  grid.rect(rep(x, ncol(z)), rep(y, each=nrow(z)), 
             width=width, height=1, default.units="native",
             gp=gpar(col=col, fill=col))
}
