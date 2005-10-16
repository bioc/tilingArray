# used for Figure 4, more than one density plot in one figure
showDens = function(z, breaks, col, ylab="", densLabels=NULL,  ...) {

  if (!is.null("densLabels")) stopifnot(length(densLabels)==length(z))
  
  y  = matrix(NA, nrow=length(breaks)-1, ncol=length(z))
  scaleFac = numeric(length(z))
  for(k in seq(along=z)) {
    h = hist(z[[k]], breaks=breaks, plot=FALSE)
    scaleFac[k] = 1/max(h$counts) * 0.92 
    y[,k] = h$counts * scaleFac[k]
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
    if (!is.null("densLabels")) text(max(pox),0.9*max(poy),densLabels[k],
                                   col=col[k],font=2,pos=2)
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
  return(scaleFac)
}

#test:
#testdat <- list(x1=rnorm(100), x2=rnorm(100,1,1), x3=rnorm(100,0.5,1))
#showDens(testdat, breaks=seq(-4,4,0.5),col=3:5, xlab="Random Numbers")

