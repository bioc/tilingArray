plot.segmentation = function(x, ...) {
  plot(x$J, 2*(1:length(x$J)))
  grid(nx=NA, ny=NULL)
}
