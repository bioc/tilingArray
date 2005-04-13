## this runs on pavlov with R-1.9.0

celdir = "Celfiles"
pdfile = "tilingBook-050413.txt"
compress = TRUE

library(affy)

## read the phenoData file
pd = read.phenoData(pdfile, fill=TRUE, as.is=TRUE,
         sep="\t", header=TRUE, comment.char="")
cat("Read", pdfile, "\n")

if(compress) {
  pData(pd)$File = paste(pd$File, ".gz", sep="")
  pat = ".cel.gz$"
} else {
  pat = ".cel$"
}

stopifnot(all(pd$use %in% c("no", "yes")))
fpd = pd[ pd$use=="yes", ]
cat("Table contains", nrow(pData(pd)), "rows, using",
    nrow(pData(fpd)), "\n")

## get the filenames in the directory
dirFiles = dir(celdir, pattern=pat)
stopifnot(all(fpd$File %in% dirFiles))

a=read.affybatch(filenames=file.path(celdir, fpd$File), 
                 phenoData=fpd, verbose=TRUE, compress=compress)

x = new("exprSet", exprs=intensity(a), phenoData=phenoData(a))
save(x, file="x.Rdata", compress=TRUE)

