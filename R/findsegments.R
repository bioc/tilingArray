findsegments <- function(x, Km, maxk) {
  G = Gmean(x, maxk)
  .Call("findsegments", G, as.integer(Km), PACKAGE="tilingArray")
}


