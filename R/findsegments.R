findSegments = function(x, maxcp, maxk, verbose=TRUE)
{
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
  if(verbose)
    cat(sprintf("findsegments: Calculation of cost matrix Gmean, n=%d, maxk=%d.\n",
                n, as.integer(maxk)))
  G = costMatrix(x, maxk)
  if (verbose) cat("Segmentation by dynamic programming.\n")
  res = .Call("findsegments", G, maxcp, verbose, PACKAGE="tilingArray")

  res$dat         <- x
  res$residuals   <- NULL
  res$chosenSegNo <- NULL
  res$confInt     <- NULL
  res$call        <- match.call()

  class(res) = c("segmentation", class(res))
  return(res)
  
}## findSegments

