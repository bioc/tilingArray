sampleStep = function(x, step) {
  if(is.integer(x))
    x=as.numeric(x)
  .Call("sampleStep", x, step, PACKAGE="tilingArray")
}


