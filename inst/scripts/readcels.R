celdir = "Celfiles"
pdfile = "tilingBook-050417.txt"

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

a = ReadAffy(filenames=fpd$File, celfile.path=celdir,  
  phenoData=fpd, verbose=TRUE, compress=TRUE)

a = new("exprSet", exprs=intensity(a), phenoData=phenoData(a))
colnames(exprs(a)) = a$File

save(a, file="a.rda", compress=TRUE)

## Normalize
jref = which(a$NucleicAcid == "DNA")
stopifnot(length(jref)==3)

normfac = rowMeans(log(exprs(a)[, jref, drop=FALSE], 2))

x = new("exprSet",
  exprs = log(exprs(a)[,-jref, drop=FALSE], 2) - normfac,
  phenoData = phenoData(a)[-jref,])

save(x, file="x.rda", compress=TRUE)

