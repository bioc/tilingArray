##--------------------------------------------------
## Create a segmentation object
##--------------------------------------------------
segment = function(y, maxseg, maxk) {
  
  if(!is.matrix(y))
    y = matrix(y, ncol=1)
  n = nrow(y)
  maxseg = as.integer(maxseg)
  maxk  = as.integer(maxk)
  if(!((length(maxseg)==1) && !is.na(maxseg) && (maxseg<=n)))
    stop(sprintf("maxseg must be an integer of length 1 between 1 and nrow(y)=%d", n))
  if(!((length(maxk)==1) && !is.na(maxk) && (maxk<=n)))
    stop(sprintf("maxk must be an integer of length 1 between 1 and nrow(y)=%d", n))
  verbose = as.integer(0)
  
  G = costMatrix(y, maxk)
  fs = .Call("findsegments", G, maxseg, verbose, PACKAGE="tilingArray")

  bp = vector(mode="list", length=nrow(fs$th))
  for(i in seq(along=bp)) {
    stopifnot( fs$th[i, i] == n+1 )
    bp[[i]] = cbind(estimate=fs$th[i, 0:(i-1)])
  }
  
  new("segmentation",
    y = y,
    breakpoints = bp,
    negloglik = fs$J,
    hasConfint = rep(FALSE, length(bp)))
} ## segment


##--------------------------------------------------
## Simple plot method
##--------------------------------------------------
setMethod("plot", "segmentation", 
  function(x, y, xlim, xlab="x", ylab="y", bpcol, ...) {
    
  validObject(x)

  y = as.integer(y)
  if(!( (length(y)==1) && (!is.na(y)) && (y>=1) && (y<=length(x@breakpoints))))
    stop(sprintf("'y' must be an integer of length 1 with values between 1 and %d.",
                 length(x@breakpoints)))

  breakp = x@breakpoints[[y]]
  ply = x@y
  plx = row(ply)
  
  ## for handling too large segmentation objects:
  if(!missing(xlim)) {
    stopifnot(is,numeric(xlim), length(xlim)==2)
    ply = ply[xlim[1]:xlim[2], ]
    plx = plx[xlim[1]:xlim[2], ]
    breakp = breakp[ (breakp[, "estimate"]>=xlim[1]) & (breakp[, "estimate"]<=xlim[2]), ]
  }

  
  if(missing(bpcol)){
    ##bpcol=hex(polarLAB(70, 35, seq(0, 360, length=y)[-1]))
    bpcol=rainbow(y)
  } else {
    if(length(bpcol)!=y)
      stop(sprintf("length of 'bpcol' is %d but must be the same as y=%d",
                   length(bpcol), as.integer(y))) 
  }

  plot(plx, ply, xlab=xlab, ylab=ylab, ...)
  abline(v=breakp[, "estimate"]-0.5, col=bpcol, lty=2) ## draw segment borders

  if (x@hasConfint[y]) { ##  confidence intervals
    for(j in 1:2) {
      at = breakp[, c("lower", "upper")[j]]
      wh = !is.na(at)
      mtext(side=1, at=at[wh]-0.5, adj=0.5, text=c("(", ")")[j], line=-0.7, col=bpcol[wh])
    }
  }
}) # plot.segmentation
    
##--------------------------------------------------
## Show method
##--------------------------------------------------
setMethod("show", "segmentation", 
  function(object) {
    ans = c(
      sprintf("Object of class 'segmentation':\n"),
      sprintf("Data matrix: %d x %d\n", nrow(object@y), ncol(object@y)),
      sprintf("Change point estimates for S = 1:%d\n", length(object@breakpoints)))
    if(any(object@hasConfint)) {
      wh = which(object@hasConfint)
      ans = c(ans,
      sprintf("Confidence intervals for %d fits from S = %d to %d\n",
              length(wh), wh[1], wh[length(wh)]))
    }
    cat(ans, "\n")
  })
