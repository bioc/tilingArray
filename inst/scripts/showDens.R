showDens = function(z, breaks, col, ylab="", ...) {
  y  = matrix(NA, nrow=length(breaks)-1, ncol=length(z))
  for(k in seq(along=z)) {
    h = hist(z[[k]], breaks=breaks, plot=FALSE)
    y[,k] = h$density/max(h$density) * 0.92
    if(k==1) {
      mids = h$mids
    } else {
      stopifnot(identical(mids, h$mids))
    }
  }

  plot(breaks[c(1, length(breaks))], c(0, length(z)), ylab=ylab, type="n", yaxt="n", ...)
  for(k in seq(along=z)) {
    poy = (length(z)-k) + c(0, rep(y[,k], each=2), 0, 0)
    pox = breaks[c(rep(seq(along=breaks), each=2), length(breaks))]
    polygon(pox, poy, col=col[k])
  }

  if(FALSE) {
    dz = lapply(z, density, from=from, to=to)
    y = sapply(dz, "[[", "y")
    x = sapply(dz, "[[", "x")
    stopifnot(all(apply(x, 1, function(p) length(unique(p)))==1))
    for(i in 1:ncol(y))
      y[,i] = y[,i]/max(y[,i]) + (i-1)
    matplot(x[,1], y, type="l", yaxt="n", lwd=3, lty=1, ...)
  }

}
