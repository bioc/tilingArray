readCel2eSet = function(filename, adf, path=".", rotated=FALSE, ...) {

  if(!missing(filename)){
    if(!missing(adf))
      stop("If 'filename' is specified, please do not specify 'adf'")
    adf = new("AnnotatedDataFrame", data=data.frame(filename=I(filename)),
        varMetadata=data.frame(labelDescription=I(c(filename="Infered from 'filename' argument of readCel2eSet"))))
  } else {
    if(missing(adf))
      stop("Please specify either 'adf' or 'filename'")
    if(!("filename" %in% varLabels(adf)))
      stop("Please let 'adf' contain a column 'filename'")
    filename = adf$filename
  }

  
  a = ReadAffy(filenames=filename, celfile.path=path, verbose=TRUE)
  ex = intensity(a)
  ex = matrix(as.integer(ex), nrow=nrow(ex), ncol=ncol(ex))
  
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

  new("ExpressionSet", exprs=ex, phenoData=adf, ...)
}

