findsegments <- function(x, maxcp, maxk) {
  G = Gmean(x, maxk)
  .Call("findsegments", G, as.integer(maxcp), PACKAGE="tilingArray")
}


