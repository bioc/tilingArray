findsegments = function(x, maxcp, maxk, verbose=0) {
  if(is.vector(x))
    x = matrix(x, nrow=length(x), ncol=1)

  n       = nrow(x)
  maxcp   = as.integer(maxcp)
  maxk    = as.integer(maxk)
  verbose = as.integer(verbose)
  if(maxcp>n)
    stop(sprintf("maxcp=%d must not be larger than nrow(x)=%d", maxcp, n))
  if(maxk>n)
    stop(sprintf("maxk=%d must not be larger than length(x)=%d", maxk, n))
  if(verbose>=2)
    cat(sprintf("findsegments: calculating Gmean, n=%d, maxk=%d.\n",
                n, as.integer(maxk)))

  G = Gmean(x[,1], maxk)
  if(ncol(x)>=2) {
    for(i in 2:ncol(x))
      G = G + Gmean(x[,i], maxk)
  }
  
  res = .Call("findsegments", G, maxcp, verbose, PACKAGE="tilingArray")
  class(res) = c("segmentation", class(res)) 
  return(res)
}


