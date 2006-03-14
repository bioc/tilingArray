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
     
  ## check that nrow(y)==length(x)==length(flag)
  n = nrow(object@y)
  if(!(length(object@x)%in%c(0,n)))
    return(FALSE)
  if(!(length(object@flag)%in%c(0,n)))
    return(FALSE)

  if(!is.na(object@nrSegments))
    if(object@nrSegments<1 && object@nrSegments>length(object@breakpoints))
      return(FALSE)
  
  isGood = TRUE
  ## check the elements of the breakpoints slot (which is a list)
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
      flag = "integer",
      breakpoints = "list",
      negloglik = "numeric",
      hasConfint = "logical",
      nrSegments = "integer"
   ),
   prototype = list(
      y = matrix(0, nrow=0, ncol=0),
      x = numeric(0),
      flag = integer(0),
      breakpoints = list(),
      negloglik = numeric(0),
      hasConfint = logical(0),
      nrSegments = as.integer(NA)
   ),
   validity = validSegmentation) ## see above
