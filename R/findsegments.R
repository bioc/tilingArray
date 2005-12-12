findSegments = function(x, maxcp, maxk, verbose=TRUE)
{
  if (verbose) cat("Assessing arguments...\n")
  if(is.matrix(x)) {
    n = nrow(x)
  } else {
    n = length(x)
  }
  maxcp   = as.integer(maxcp)
  maxk    = as.integer(maxk)
  verbose = as.integer(verbose)
  if(maxcp>n)
    stop(sprintf("maxcp=%d must not be larger than nrow(x)=%d", maxcp, n))
  if(maxk>n)
    stop(sprintf("maxk=%d must not be larger than length(x)=%d", maxk, n))
  if(verbose>=2)
    cat(sprintf("findsegments: calculating Gmean, n=%d, maxk=%d.\n",
                n, as.integer(maxk)))
  if (verbose) cat("Computing cost matrix for segmentation...\n")
  G = costMatrix(x, maxk)
  if (verbose) cat("Running Picard's segmentation algorithm...\n")
  res = .Call("findsegments", G, maxcp, verbose, PACKAGE="tilingArray")

  res$dat         <- x
  res$residuals   <- NULL
  res$chosenSegNo <- NULL
  res$confInt     <- NULL
  res$call        <- match.call()

  class(res) = c("segmentation", class(res))
  return(res)
  
}#findSegments


plot.segmentation = function(x, nSegments, bcol=NULL, from=NULL, to=NULL, ...){

  stopifnot(all(c("th","dat") %in% names(x)), is.numeric(x$dat))
  
  y   <- as.vector(x$dat)
  Index <- rep(1:nrow(as.matrix(x$dat)), length.out=length(y))
  
  if (!is.null(nSegments)){
    stopifnot(is.numeric(nSegments), nSegments>1, nSegments<= nrow(x$th))
  } else { # WORKING SOLUTION: use largest drop of RSS to determine nSeg
    if (is.null(x$residuals))
      stop("\nObject does not contain computed residuals!\nUse function 'confint' to compute those or specify\nnumber of segments in segmentation to plot.\n")
    rss     <- sapply(x$residuals, function(z) sum(z^2))
    if (length(rss)<2)
      nSegments <- x$chosenSegNo[1]
    else {
      rssdiff <- diff(rss)
      nSegments <- x$chosenSegNo[which.min(rssdiff)+1]
    }
  }#else

  breakp <- x$th[nSegments, 1:(nSegments-1), drop=TRUE]
  ncp    <- nSegments - 1 # number of change points
  
  if (is.null(x$confInt))
    ci <- matrix(breakp, nrow=length(breakp), ncol=3, byrow=FALSE)
  else
    ci <- x$confInt[[match(nSegments,x$chosenSegNo)]] # confidence intervals

  
  # for handling too large segmentation objects:
  if (!is.null(from)|!is.null(to)){
    if (is.null(from)) from <- 1
    if (is.null(to))   to   <- length(y)
    stopifnot(from >= 1, from < to, to <= length(y))
    y          <- y[from:to]
    Index      <- Index[from:to]
    keep.index <- which( breakp>from & breakp<=to)
    breakp     <- breakp[keep.index]
    ci         <- ci[keep.index, ,drop=FALSE]
  }

  plot(x=Index, y=y, ...)
  
  if (is.null(bcol))
    mycols <- 1:length(breakp) + 1
  else
    mycols <- bcol
  abline(v=ci[,2]-0.5,col=mycols,lty=2) # draw segment borders

  if (!is.null(x$confInt)){ # draw change point confidence intervals
    mtext(side=1, at=ci[,1]-0.5, text=rep("(",ncp), line=-0.7, col=mycols)
    mtext(side=1, at=ci[,3]+0.5, text=rep(")",ncp), line=-0.7, col=mycols)
  }
  
  invisible(list(breakp=breakp, confInt=ci))
} # plot.segmentation
    
