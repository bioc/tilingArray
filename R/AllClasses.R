## ==========================================================================
## (C) W. Huber 2005
## segmentation: a class to contain segmentation models (i.e. fits of
## piecewise constant functions)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##--------------------------------------------------
## Validity method
##--------------------------------------------------
validSegmentation = function(object) {

  ## check that slots 'breakpoints', 'negloglik', 'hasConfint' all have the same length
  if(!((length(object@breakpoints)==length(object@negloglik)) &&
       (length(object@breakpoints)==length(object@hasConfint))))
    return(FALSE)
     
  ## check that nrow(y)==length(x)
  n = nrow(object@y)
  if(!is.null(object@x) && (length(object@x)!=n))
    return(FALSE)
  
  isGood = TRUE
  for(i in seq(along=object@breakpoints)) {
    b = object@breakpoints[[i]]
    if(!((nrow(b)==i-1) && (ncol(b) %in% c(1,3)) &&
         ("estimate" %in% colnames(b)) && all(b[, "estimate"] <= n, na.rm=TRUE)))
      isGood = FALSE
  }
  return(isGood)
}

##--------------------------------------------------
## Definition
##--------------------------------------------------
setClass("segmentation",
   representation(
      y = "matrix",
      x = "numeric",             
      breakpoints = "list",
      negloglik = "numeric",
      hasConfint = "logical"
   ),
   prototype = list(
      y = matrix(0, nrow=0, ncol=0),
      x = numeric(0),
      breakpoints = list(),
      negloglik = numeric(0),
      hasConfint = logical(0)
   ),
   validity = validSegmentation) ## see above
