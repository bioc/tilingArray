##-----------------------------------------------------------------------
## This function calculates the cost matrix for a homoscedastic model
## G[k,i] is the sum of squared residuals of a segment from i to (i+k-1)
##-----------------------------------------------------------------------
costMatrix = function(x, maxk) {
  if(!is.numeric(maxk)||maxk<=1||length(maxk)!=1)
    stop("'maxk' must be a single positive number.")
  if(!is.numeric(x)||!(is.vector(x)||is.matrix(x)))
    stop("'x' must be a numeric vector or matrix.")

  if(is.vector(x))  {
    r = x
    q = x*x
    S = 1
  } else {
    r = rowSums(x)
    q = rowSums(x*x)
    S = ncol(x)
  }

  ## see Wolfgang's handwritten notes for explanation of the algebra
  c = cumsum(r)
  d = cumsum(q)
  
  G = matrix(as.numeric(NA), nrow=maxk, ncol=n)
  for(k in 1:maxk) {
    i   = 0:(n-k)
    dki = d[i+k]-d[i]
    cki = c[i+k]-c[i]
    G[k, ] = (dki - cki*cki/S)/S
  }
  return(G)
}

