##-----------------------------------------------------------------------
## This function calculates the cost matrix for a homoscedastic model
## G[k,i] is the sum of squared residuals of a segment from i to (i+k-1)
##-----------------------------------------------------------------------
costMatrix = function(x, maxk) {
  if(!is.numeric(x)||!(is.vector(x)||is.matrix(x)))
    stop("'x' must be a numeric vector or matrix.")

  if(is.vector(x))  {
    r = x
    q = x*x
    d = 1
  } else {
    r = rowSums(x)
    q = rowSums(x*x)
    d = ncol(x)
  }
  n = length(r)
  
  if(!is.numeric(maxk)||(length(maxk)!=1)||(maxk<1)||(maxk>n))
    stop(sprintf("'maxk' must be a single number between 1 and the number of rows of 'x': %d.", n))

  ## see inst/doc/costMatrix.tex (.pdf) for explanation of the algebra
  cr = cumsum(r)
  cq = cumsum(q)
  
  G = matrix(as.numeric(NA), nrow=maxk, ncol=n)
  k = 1:maxk
  G[, 1] = cq[k] - cr[k]*cr[k]/(k*d)
  for(k in 1:(maxk-1)) {
    i   = 1:(n-k)
    j   = 2:(n-k+1)
    cqk = cq[i+k]-cq[i]
    crk = cr[i+k]-cr[i]
    G[k,j] = cqk - crk*crk/(k*d)
  }
  return(G)
}

