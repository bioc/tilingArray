otherStrand = function(x) {
  x  = as.character(x)
  mt = match(x, c("-", "+"))
  if(any(is.na(mt)))
    stop("All elements of 'x' must be + or -.")
  return(c("+", "-")[mt])
}
