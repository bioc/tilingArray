library("Biobase")
library("RColorBrewer")

if(!exists("a")) {
  cat("Loading a.rda\n")
  load("a.rda")
}
if(!exists("probeAnno")) {
  cat("Loading probeAnno.rda\n")
  load("probeAnno.rda")
}
if(!exists("xn")) {
  cat("Loading xn.rda\n")
  ## load("seg-polyA-050525/xn.rda")
  load("seg-tot-050525/xn.rda")
}

if(!exists("refSigThresh")) {
  chrstr = paste(rep(1:17, each=2),
    rep(c("+", "-"), 17), sep=".")
  allPM = unique(unlist(lapply(chrstr, function(chr)
    get(paste(chr, "index", sep="."), probeAnno))))
  refSigThresh = quantile(refSig[allPM], probs=0.05)
}
cat("refSigThresh=", signif(refSigThresh, 3), "\n")

jref = which(a$NucleicAcid == "DNA")
stopifnot(length(jref)==3)

fn = c("05_04_27_2xpolyA_NAP3.cel.gz",
      "05_04_26_2xpolyA_NAP2.cel.gz",
      "05_04_20_2xpolyA_NAP_2to1.cel.gz")

## UFD2: Chr 4-, 118-122000
##sta = probeAnno$"4.-.start"
##ind = probeAnno$"4.-.index"
##sel = (sta>=115000 & sta<=125000)

sta = probeAnno$"1.-.start"
ind = probeAnno$"1.-.index"
sel = (sta>=5000 & sta<=11000)

isel = ind[sel]
xsel = sta[sel]

y1 = rowMeans(log(exprs(a)[isel, jref], 2))
y2 = rowMeans(log(exprs(a)[isel, fn], 2))
y3 = rowMeans(exprs(xn)[isel,  ])
y4 = y3
y4[refSig[isel]<refSigThresh] = NA
ysel = cbind(y1, y2, y2-y1, y3, y4)

colnames(ysel) = c("DNA", "unnormalized", "normalization method 1",
      "normalization method 2","normalization method 3")

graphics.off()
x11(width=8, height=12)
par(mfrow=c(ncol(ysel), 1))
rgnf = range(refSig)
cols = colorRamp(c("red", "yellow", "blue"))((refSig[isel]-rgnf[1])/(rgnf[2]-rgnf[1])) / 256
## cols = colorRamp(brewer.pal(11, "Spectral"))((refSig[isel]-rgnf[1])/(rgnf[2]-rgnf[1])) / 256
cols = rgb(cols[,1], cols[,2], cols[,3])
## cols = rainbow(256)[ ceiling(rank(refSig[isel])/length(isel) * 256) ]

for(i in 1:ncol(ysel)) {
  py  = ysel[,i]
  sl  = !is.na(py)
  plot(xsel[sl], py[sl], col=cols[sl], pch=20, xlab="coordinates", ylab="signal",
       main=colnames(ysel)[i])
}

dev.copy(pdf, file="assessNorm.pdf", width=8, height=12); dev.off()

## x11(width=6, height=6)
## plot(ysel[, c("RNA new norm.", "RNA old norm.")], pch=17, col=c("black", "red")[1+badProbes[isel]])


stop()

criterion = function(x) {
  x  = x[!is.na(x)]
  dd = abs(x[-(1:4)]-x[(0:3)-length(x)]) / IQR(x, na.rm=TRUE)
  quantile(dd, probs=c(0.5))
}

sel = (sta>=100000 & sta<=140000)
isel = ind[sel]
xsel = sta[sel]
y1 = rowMeans(log(exprs(a)[isel, jref], 2))
y2 = rowMeans(log(exprs(a)[isel, fn], 2))
y3 = rowMeans(exprs(xn)[isel,  ])
y4 = y3
y4[refSig[isel]<refSigThresh] = NA
ysel = cbind(y1, y2, y2-y1, y3, y4)
print(apply(ysel[order(xsel), ], 2, criterion))

 
