##-----------------------------------------------------------------------
## This function calculates the cost matrix for a homoscedastic model
## G[k,i] is the sum of squared residuals of a segment from i to (i+k-1)
##-----------------------------------------------------------------------
Gmean = function(x, maxk) {
  if(!is.numeric(maxk)||maxk<=1||length(maxk)!=1)
    stop("'maxk' must be a single positive number.")
  if(!is.numeric(x)||!is.vector(x))
    stop("'x' must be a numeric vector.")
  
  n  = length(x)
  cx = cumsum(x)
  x2 = x*x
  A = G = matrix(as.numeric(NA), nrow=maxk, ncol=n)
  A[1, ] = x
  A[, 1] = cx[1:maxk]
  i = 2:n
  for(k in 2:maxk) {
    A[k, i] = cx[i-1+k] - cx[i-1]
  }
  A = A*A
  G[1, ] = 0
  i = 1:n
  for(k in 2:maxk) {
    G[k,] = x2[i-1+k] - A[k,]/k + A[k-1,]/(k-1) + G[k-1,]
  }
  return(G)
}

