grid.image = function(x, y, z, uniq=uniq,
  colRamp=colorRamp(brewer.pal(9, "YlGnBu")),
  width=22, just=c("left","centre")) {
  
  mcol = (colRamp(z))/256
  mcol[uniq!=0,] = 1
  indna <- is.na(mcol[,1])
  mcol[indna,] = 1

  col = rgb(mcol[,1], mcol[,2], mcol[,3])
  grid.rect(rep(x, ncol(z)), rep(y, each=nrow(z)),just=just,
            width=width, height=1, default.units="native",
            gp=gpar(col=col, fill=col))
}
# the c("left","center") justification is important. The first one
#  is the x-axis (chromosomal coordinate) one.
#  Thus "left" is correct, if the given x-coordinates are the START positions
#  of probes and the segmenation lines should be to the left of them 
