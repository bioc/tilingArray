findsegments <- function(x, maxcp, maxk) {
  G = Gmean(x, maxk)
  res = .Call("findsegments", G, as.integer(maxcp), PACKAGE="tilingArray")
  class(res) = c("segmentation", class(res)) 
  return(res)
}


