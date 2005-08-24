confint.segmentation <- function(object,  parm="breakpoints", level = 0.95,
                                 nSegments=NULL, 
                                 het.reg = FALSE,het.err = FALSE, ...)
{
  # function to compute confidence intervals after a segmentation for a
  #  specified or all considered number of segments
  require(strucchange)

  # check arguments:
  stopifnot(#class(object)=="segmentation",
            all(c("dat","th") %in% names(object)), nrow(object$th)>1,
            is.null(nSegments)| (is.numeric(nSegments)& all(nSegments>1) &
                                 all(nSegments<=nrow(object$th))),
            is.numeric(level), level>0, level<1)
  
  x <- as.matrix(object$dat)
  n <- nrow(x)
  lx <- length(x)

  if (is.null(nSegments)) # number of Segments not specified -> use all in 'th'
    nSegments <- 2:nrow(object$th)

  allConfInt <- list()
  allResid   <- list()

  # compute confidence intervals for all specified segmentations:
  for (thisnseg in nSegments){
    
    breakp <- object$th[thisnseg, 1:(thisnseg-1), drop=TRUE]
    bp <- c(0,breakp,(n+1))
  
    segcut <- cut(1:n, bp, right=FALSE)
      
    resid <- apply(x, 2, function(z)
        unlist(tapply(z, segcut, function(this) this-mean(this)),use.names=FALSE))
    # cut the samples into segments, compute mean ms for each segment and
    #  subtract ms from each value in the segment, do this separately for each
    #  sample (column of dat) to account for sample-specific segment means
    stopifnot(length(resid)==lx)

    xobject <- list(y=x, breakpoints=breakp,residuals=resid)

    xconfint <- computeConfInt(xobject, level = level,
                               het.reg= het.reg, het.err= het.err)
    
    allConfInt <- c(allConfInt, list(xconfint))
    allResid   <- c(allResid, list(resid))
    
  }# for (thisnseg in nSegments)
  
  #return(xconfint)
  # add results to segmentation object and return:
  object$chosenSegNo   <- nSegments
  names(allConfInt) <- names(allResid) <- paste(nSegments,"Segments",sep=" ")
  object$confInt       <- allConfInt
  object$residuals     <- allResid
  object$call          <- c(object$call, match.call())
  return(object)
}
