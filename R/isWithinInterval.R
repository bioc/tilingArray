isWithinInterval = function(x, start, end) {
  if(!(is.numeric(x) && is.vector(x)))
    stop("'x' must be a numeric vector.")
  if(!(is.numeric(start) && is.vector(start) && is.numeric(end) && is.vector(end) &&
       length(end)==length(start) && identical(names(start), names(end))))
    stop("'start' and 'end' must be numeric vectors of the same length and with same names.")
  n = length(x)
  m = length(start)

  mx  = matrix(x, nrow=n, ncol=m)
  res = ((mx >= matrix(start, nrow=n, ncol=m, byrow=TRUE)) &
         (mx <= matrix(end,   nrow=n, ncol=m, byrow=TRUE)))

  rownames(res)=names(x)
  colnames(res)=names(start)
  return(res)
}
