readCel2eSet = function(filenames, path=".", rotated=FALSE) {
  
  a = ReadAffy(filenames=filenames, celfile.path=path, verbose=TRUE)
  ex = intensity(a)
  
  if(!rotated) {
    n = a@nrow
    if(a@ncol!=n)
      stop("Don't know how to deal with chips for which ncol != nrow")
    
    xy2i = function(coord) 1 + coord[,1] + n*coord[,2]
    i2xy = function(i)  cbind( x = (i-1) %% n, y = (i-1) %/% n )
    rotate = function(coord) {
      off = (n-1)/2
      rot = cbind(c(0,1), c(-1,0))
      off + (coord-off) %*% rot 
    }
    ex = ex[ xy2i(rotate(i2xy(1:(n*n)))),,drop=FALSE ]
  }
  
  e = new.env(parent = baseenv())
  assign("exprs", ex, e)
  
  new("eSet", assayData=e, phenoData=phenoData(a), sampleNames=sampleNames(a))
}

