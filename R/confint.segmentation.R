## (c) wolfgang huber 2005
## compute confidence intervals of segmentation(s)
setMethod("confint", "segmentation",
  function(object, parm, level=0.95, het.reg = FALSE, het.err = FALSE, ...)
{
  validObject(object)
  ## the list with breakpoints
  bpL = object@breakpoints

  if(length(bpL)<2)
    stop("length of breakpoints list must be >=2")
  
  ## Check arguments
  if(!missing(parm)) {
    parm = as.integer(parm)
    if(!(is.numeric(parm) &&
         all(parm>1) &&
         all(parm<=length(bpL))))
      stop(sprintf("'parm' must be numeric with values between 2 and %d", length(bpL)))
  } else {
    parm = 2:length(bpL)
  }

  ## transpose: this way the replicate data points (different columns of object@y)
  ##   come after another.
  toy = t(object@y)
    
  ## loop over parm
  for (j in parm) {
    ## Breakpoints:
    ## -- subtract one, since "object@breakpoints" stores the
    ## indices of last points of each segment plus 1 (i.e. the first point of
    ## the next segment), while strucchange expects indices of last points.
    ## -- multiply by nrow(toy) since this is the number of replicates at each
    ## x-position.
    bpj = bpL[[j]][, "estimate"]
    breaks = (bpj-1) * nrow(toy)

    ## Residuals (this really ought to be calculated by the fit function,
    ##  but currently it isn't)
    res = toy
    bp = c(0, breaks, length(toy))
    for(k in 2:length(bp)) {
      rg = (bp[k-1]+1):bp[k]
      res[rg] = res[rg]-mean(res[rg])
    }
    browser()
    
    ## Assemble a pretend "breakpointsfull" object, with which we can call 
    ## "confint.breakpointsfull" from the strucchange package
    bpp = list(X = matrix(as.integer(1), nrow=length(toy), ncol=1),
               y = toy,
               nobs = length(toy),
               breaks = breaks,
               res = res,
               nreg = NULL, datatsp = NULL)
    class(bpp) = "breakpointsPretend"
    
    ci =  strucchange:::confint.breakpointsfull(bpp, level=level, het.reg=het.reg, het.err=het.err, ...)

    ## extract the confidence intervals and add back 1 to be consistent with our ways
    m = ci$confint
    stopifnot(!is.null(m), all(m[, 2]==breaks))
    m[] = as.integer(round(m/nrow(toy) + 1))
    colnames(m) = c("lower", "estimate", "upper")
    bpL[[j]] = m
    
  } ## for j
  object@breakpoints = bpL
  object@hasConfint[parm] = TRUE
  if(length(parm)==1)
    object@nrSegments=parm
  
  return(object)
})

## argument "breaks" is ignored, it is NULL anyway
residuals.breakpointsPretend = function(object, breaks, ...) {
  stopifnot(is.null(breaks))
  object$res
}

breakpoints.breakpointsPretend = function(obj, breaks, ...) {
  stopifnot(is.null(breaks))
  list(breakpoints =obj$breaks)
}
