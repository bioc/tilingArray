## compute confidence intervals of segmentation(s)
confint.segmentationNew <- function(object,  parm="breakpoints", level = 0.95,
                                 nSegments, 
                                 het.reg = FALSE, het.err = FALSE, ...)
{
  ## Check arguments
  ##
  ## FIXME: when we have an S4 implementation of segmentation objects,
  ##   this kind of checking can and should be done in the class'
  ##   validObject function.
  if(!(inherits(object, "segmentation") &&
       all(c("dat", "th") %in% names(object)) &&
       is.matrix(object$th) &&
       is.matrix(object$dat) &&
       (nrow(object$th)>1)) &&
       (ncol(object$dat)>=1) &&
       (nrow(object$dat)==object$th[length(object$th)]-1))
    stop("'object' is not a valid instance of class 'segmentation'")
  
  if(!missing(nSegments))
    if(!(is.numeric(nSegments) &&
         all(nSegments>1) &&
         all(nSegments<=nrow(object$th))))
      stop("'nSegments' must be numeric with values between 2 and nrow(object$th)")

  if(!(is.numeric(level) && (length(level)==1) && abs(level-0.5)<0.5))
      stop("'nSegments' must be numeric of length 1 with values between 0 and 1")
  
  ## transpose: this way the replicate data points (different columns of object$dat)
  ##   come after another.
  y = t(object$dat)
  
  ## if number of Segments not specified, use all in 'th'
  if (missing(nSegments)) 
    nSegments <- 2:nrow(object$th)

  allConfInt = vector(mode="list", length=length(nSegments))
  names(allConfInt) = paste(nSegments)
  
  ## loop over nSegments
  for (j in seq(along=nSegments)) {

    ## Breakpoints: subtract one, since "object$th" stores the indices of
    ## last points of each segment plus 1 (i.e. the first point of the
    ## next segment), while strucchange expects indices of last points.
    ## Multiply by nrow(y) since this is the number of replicates at each
    ## x-position.
    ## Omit the last breakpoint since that is implicit in strucchange.
    breaks = (object$th[nSegments[j], 1:(nSegments[j]-1)] - 1) * nrow(y)

    ## Residuals
    res = y
    bp = c(0, breaks, length(y))
    for(k in 2:length(bp)) {
      rg = (bp[k-1]+1):bp[k]
      res[rg] = res[rg]-mean(res[rg])
    }

    ## Assemble a pretend "breakpointsfull" object, with which we can call 
    ## "confint.breakpointsfull" from the strucchange package
    bpp = list(X=matrix(as.integer(1), nrow=length(y), ncol=1),
               y=y,
               nobs=length(y),
               breaks=breaks,
               res=res,
               nreg=NULL, datatsp=NULL)
    class(bpp) = "breakpointsPretend"
    
    ci = confint.breakpointsfull(bpp, level = level, het.reg= het.reg, het.err= het.err)

    ## extract the confidence intervals and add back 1 to be consistent with our ways
    stopifnot("confint" %in% names(ci), all(ci$confint[,2]==breaks))
    allConfInt[[j]] = ci$confint+1
  } ## for j
  
  ## add results to segmentation object and return:
  object$nSegments = nSegments   ## FIXME: this is redundant with the names of the confint slot - omit?
  object$confint   = allConfInt
  object$call      = c(object$call, match.call())
  return(object)
}

## argument "breaks" is ignored, it is NULL anyway
residuals.breakpointsPretend = function(x, breaks) {
  stopifnot(is.null(breaks))
  x$res
}
breakpoints.breakpointsPretend = function(x, breaks) {
  stopifnot(is.null(breaks))
  list(breakpoints = x$breaks)
}
