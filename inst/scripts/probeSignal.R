library("Biobase")
library("prada")
source("colorRamp.R")

if(!exists("x"))load("x.Rdata")
if(!exists("probeAnno"))load("probeAnno.rda")

iPM = lapply(paste(1:16), function(chr) {
  lapply(c("+", "-"), function(s) {
    get(paste(chr, s, "index", sep="."), probeAnno)
  })
})

iPM = unlist(iPM)

graphics.off()
par(mfcol=c(2,3))

DNAhybes = c("09_11_04_S96_genDNA_16hrs_45C_noDMSO.cel.gz",
  "041119_S96genDNA_re-hybe.cel.gz",           
  "041120_S96genDNA_re-hybe.cel.gz")

for (fn in DNAhybes) {
  j = match(fn, x$File)
  stopifnot(!is.na(j))
  hist(log(exprs(x)[iPM, j], 2), main=paste("PM", fn), breaks=50, col="lightblue")
  hist(log(exprs(x)[iPM+2560, j], 2), main=paste("MM", fn), breaks=50, col="mistyrose")
}

dev.copy(pdf, file="probeSignal-hist.pdf", width=12, height=6)
dev.off()


par(mfcol=c(1,3))
for (j in 2:length(DNAhybes)) {
  for (i in 1:(j-1)) {
    smoothScatter(log(exprs(x)[iPM, c(i,j)], 2))
  }
}
dev.copy(pdf, file="probeSignal-scp.pdf", width=12, height=4.6)
dev.off()


