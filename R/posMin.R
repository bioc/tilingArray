posMin = function(x, ...) {
  x=x[x>=0]
  if(length(x)>0) {
    min(x, ...)
  } else {
    as.numeric(NA)
  }
}

