findsegments <- function(x, maxcp, maxk, verbose=0) {
  if(verbose>=2)
    cat(sprintf("findsegments: calculating Gmean, n=%d, maxk=%d.\n",
                length(x), as.integer(maxk))) 
  G = Gmean(x, maxk)
  res = .Call("findsegments", G, as.integer(maxcp), as.integer(verbose), PACKAGE="tilingArray")
  class(res) = c("segmentation", class(res)) 
  return(res)
}


