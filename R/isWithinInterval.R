isWithinInterval = function(x, start, end) {
  n = length(x)
  m = length(start)

  if(!(is.numeric(x) && is.vector(x) && n>0))
    stop("'x' must be a numeric vector of length > 0.")
  
  if(!(is.numeric(start) && is.vector(start) && is.numeric(end) &&
       is.vector(end) && length(end)==m &&
       identical(names(start), names(end)) && m>0))
    stop(paste("'start' and 'end' must be numeric vectors of the",
         "same length > 0, and with same names."))

  if(!require("SparseM")) {
    warning(paste("Package 'SparseM' could not be loaded, reverting",
        "to normal matrices (this maybe rather memory inefficient)."))
   
    mx  = matrix(x, nrow=n, ncol=m)
    res = ((mx >= matrix(start, nrow=n, ncol=m, byrow=TRUE)) &
           (mx <= matrix(end,   nrow=n, ncol=m, byrow=TRUE)))
  } else {

    
  }
  
  rownames(res)=names(x)
  colnames(res)=names(start)
  return(res)
}
