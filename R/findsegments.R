findSegments = function(x, maxcp, maxk, verbose=0) {
  if(is.matrix(x)) {
    n = nrow(x)
  } else {
    n = length(x)
  }
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

  G = costMatrix(x, maxk)
  res = .Call("findsegments", G, maxcp, verbose, PACKAGE="tilingArray")
  class(res) = c("segmentation", class(res)) 
  return(res)
}


