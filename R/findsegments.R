findsegments = function(x, maxcp, maxk, verbose=0) {
  maxcp   = as.integer(maxcp)
  maxk    = as.integer(maxk)
  verbose = as.integer(verbose)
  if(maxcp>length(x))
    stop(sprintf("maxcp=%d must not be larger than length(x)=%d", maxcp, length(x)))
  if(maxk>length(x))
    stop(sprintf("maxk=%d must not be larger than length(x)=%d", maxk, length(x)))
  if(verbose>=2)
    cat(sprintf("findsegments: calculating Gmean, n=%d, maxk=%d.\n",
                length(x), as.integer(maxk)))
  
  G = Gmean(x, maxk)
  gc()
  res = .Call("findsegments", G, maxcp, verbose, PACKAGE="tilingArray")
  class(res) = c("segmentation", class(res)) 
  return(res)
}


