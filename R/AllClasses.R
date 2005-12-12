## ==========================================================================
## (C) W. Huber 2005
## segmentation: a class to contain segmentation models (i.e. fits of
## piecewise constant functions)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##--------------------------------------------------
## Validity method
##--------------------------------------------------
validSegmentation = function(object) {
  if(!((length(object@breakpoints)==length(object@negloglik)) &&
       (length(object@breakpoints)==length(object@hasConfint))))
    return(FALSE)
     
  n = nrow(object@y)
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
      y  = "matrix",
      breakpoints = "list",
      negloglik = "numeric",
      hasConfint = "logical"
   ),
   prototype = list(
      y  = matrix(as.numeric(NA), nrow=0, ncol=0),
      breakpoints = list(),
      negloglik = numeric(0),
      hasConfint = logical(0)
   ),
   validity = validSegmentation) ## defined in methods-segmentation.R
