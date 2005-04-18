## do the probes that have no or weak DNA signal cluster in particular
## regions or are they sporadic?
library("Scerevisiaetilingprobe")
## if(!exists("x"))load("x.rda")
if(!exists("a"))load("a.rda")
if(!exists("probeAnno"))load("probeAnno.rda")

if(!exists("dnasig")) {
  jref = which(a$NucleicAcid == "DNA")
  stopifnot(length(jref)==3)
  dnasig = rowMeans(log(exprs(a)[, jref, drop=FALSE], 2))
}

if(!exists("bc"))
  bc = basecontent(Scerevisiaetilingprobe$sequence)

for(chr in 1:16) {
  for(i in 1:2) {
    s = c("+", "-")[i]
    ind = get(paste(chr, s, "index", sep="."), probeAnno) 
    sta = get(paste(chr, s, "start", sep="."), probeAnno)
    y   = dnasig[ind]
    d   = density(y, adjust=0.1)
    switch(i,
        plot(d, main=chr),
        lines(d$x, d$y+0.1*diff(range(d$y)), col="blue"))
    percent = percent + sum(iw)/length(iw)*50
  }
  cat(sprintf("%2d: ", as.integer(chr)), round(percent, 1), "%", "\n", sep="")
}


## This code produces density of "eliminated" probes

isweak = is.na(exprs(x)[,1])
par(mfrow=c(4,4))
for(chr in 1:16) {
  percent = 0
  for(i in 1:2) {
    s = c("+", "-")[i]
    ind = get(paste(chr, s, "index", sep="."), probeAnno) 
    sta = get(paste(chr, s, "start", sep="."), probeAnno)
    iw  = isweak[ind]
    d   = density(sta[iw], adjust=0.1)
    switch(i,
        plot(d, main=chr),
        lines(d$x, d$y+0.1*diff(range(d$y)), col="blue"))
    percent = percent + sum(iw)/length(iw)*50
  }
  cat(sprintf("%2d: ", as.integer(chr)), round(percent, 1), "%", "\n", sep="")
}

dev.copy(pdf, file="weakProbesCluster.pdf", width=12, height=10)
dev.off()
