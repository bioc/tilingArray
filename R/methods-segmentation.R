##
## methods associated with object of class "segmentation"
##   please see also file confint.segmentation.R
##
## (C) W. Huber 2005
##

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
    x = numeric(0), flag=integer(0),
    breakpoints = bp,
    logLik = fs$J,
    hasConfint = rep(FALSE, length(bp)))
} ## segment


##--------------------------------------------------
## Simple plot method
##--------------------------------------------------
setMethod("plot", "segmentation", 
  function(x, y, xlim, xlab="x", ylab="y", bpcol="black", bplty=1, pch=16, ...) {
    
  validObject(x)

  y = as.integer(y)
  if(!( (length(y)==1) && (!is.na(y)) && (y>=1) && (y<=length(x@breakpoints))))
    stop(sprintf("'y' must be an integer of length 1 with values between 1 and %d.",
                 length(x@breakpoints)))

  ply = x@y
  if(ncol(ply)>1) {
    note("slot 'y' has more than one column, calculating 'rowMeans'")
    ply = matrix(rowMeans(ply), ncol=1)
  }
  
  plx = if(length(x@x)==0) row(ply) else x@x

  ## extract the y-th set of breakpoints, and transform from index counts to actual
  ##   regressor values (plx)
  breakp = x@breakpoints[[y]]
  for(i in 1:ncol(breakp))
    breakp[,i] = plx[ breakp[,i] ]

  ## since y is the number of segments = 1 +  number of breakpoints
  stopifnot(nrow(breakp)==y-1)
  
  ## for handling large segmentation objects:
  if(!missing(xlim)) {
    stopifnot(is.numeric(xlim), length(xlim)==2)
    sel = (plx>=xlim[1]) & (plx<=xlim[2])
    ply = ply[sel, ]
    plx = plx[sel]
    breakp = breakp[ (breakp[, "estimate"]>=xlim[1]) & (breakp[, "estimate"]<=xlim[2]),, drop=FALSE ]
  }
  
  if(missing(bpcol))
    bpcol = rep(brewer.pal(9, "Set1"), ceiling(nrow(breakp)/9))[1:nrow(breakp)]
  
  plot.default(plx, ply, xlab=xlab, ylab=ylab, pch=pch, ...)
  abline(v=breakp[, "estimate"]-0.5, col=bpcol, lty=bplty) ## draw segment boundaries

  if (x@hasConfint[y]) { ##  confidence intervals
    for(j in 1:2) {
      at = breakp[, c("lower", "upper")[j]]
      isnotna = !is.na(at)
      isdup   = duplicated(at[isnotna])
      mtext(side=1, at=at[isnotna]-0.5+0.25*isdup, adj=0.5, text=c("(", ")")[j], line=-0.7, col=bpcol[isnotna])
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
      sprintf("Change point estimates for number of segments S = 1:%d\n", length(object@breakpoints)))
    if(any(object@hasConfint)) {
      wh = which(object@hasConfint)
      ans = c(ans,
      sprintf("Confidence intervals for %d fits from S = %d to %d\n",
              length(wh), wh[1], wh[length(wh)]))
    }
    ans = c(ans, 
      sprintf("Selected S = %d\n", object@nrSegments))
    cat(ans, "\n", sep="")
  })

##--------------------------------------------------
## logLik
##--------------------------------------------------
setMethod("logLik", "segmentation", 
  function(object, penalty="none", ...) {
    ll   = object@logLik
    npar = 2*(1:length(ll))
    ndat = length(object@y)
    switch(penalty,
           none = ll,
           AIC = ll - npar,
           BIC = ll - log(ndat)*npar/2,
           stop(sprintf("Invalid value '%s' for argument 'penalty'.", penalty)))

  })

