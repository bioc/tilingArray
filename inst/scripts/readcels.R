celdir = "Celfiles"
pdfile = "tilingBook-050413.txt"

library(affy)

## Read the phenoData file
pd = read.phenoData(pdfile, fill=TRUE, as.is=TRUE,
         sep="\t", header=TRUE, comment.char="")
cat("Read", pdfile, "\n")

## CEL files are compressed:
pData(pd)$File = paste(pd$File, ".gz", sep="")
pat = ".cel.gz$"

stopifnot(all(pd$use %in% c("no", "yes")))
fpd = pd[ pd$use=="yes", ]
cat("Table contains", nrow(pData(pd)), "rows, using",
    nrow(pData(fpd)), "\n")

## Get the filenames in the directory
dirFiles = dir(celdir, pattern=pat)
stopifnot(all(fpd$File %in% dirFiles))

stop()

a = read.affybatch(filenames=file.path(celdir, fpd$File), 
                 phenoData=fpd, verbose=TRUE, compress=compress)

a = new("exprSet", exprs=intensity(a), phenoData=phenoData(a))
save(a, file="a.Rdata", compress=TRUE)

## Normalize
jref = which(x$Hybe %in% c(12,23,24))
  stopifnot(length(jref)==3, !any(is.na(k)))
  return(log(exprs(x)[,k,drop=FALSE], 2) - rowMeans(log(exprs(x)[,jref], 2)))
}


x =  
save(x, file="x.Rdata", compress=TRUE)

