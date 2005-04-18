## do the probes that have no or weak DNA signal cluster in particular
## regions or are they sporadic?
library("Scerevisiaetilingprobe")
library("locfit")

## if(!exists("x"))load("x.rda")
if(!exists("a"))load("a.rda")
if(!exists("probeAnno"))load("probeAnno.rda")

if(!exists("dnasig")) {
  jref = which(a$NucleicAcid == "DNA")
  stopifnot(length(jref)==3)
  dnasig = rowMeans(log(exprs(a)[, jref, drop=FALSE], 2))
}

if(!exists("at")) {
  basc = basecontent(Scerevisiaetilingprobe$sequence)
  at = basc[,"A"]+basc[,"T"]
}

mylf = function(x, y, i, ...) {
  lf = locfit.raw(x, y, alpha=0.02, deg=1)
  rx = range(x)
  px = seq(rx[1], rx[2], length=400)
  plf = predict(lf, newdata=px)
  list(x=px, y=plf)
}

doPDF = TRUE
if(doPDF){
  pdf(file="weakProbesCluster.pdf", width=10, height=12)
  par(mfcol=c(6,2), mai=c(0.5, 0.5, 0.2, 0.1))
  nrchr=16
} else {
  graphics.off(); x11(width=12, height=10)
  par(mfcol=c(4,2))
  nrchr=4
}

for(chr in 1:nrchr) {
  for(i in 1:2) {
    s = c("+", "-")[i]
    ind  = get(paste(chr, s, "index", sep="."), probeAnno) 
    sta  = get(paste(chr, s, "start", sep="."), probeAnno)
    cat(chr, s, "  ", sep="")
    
    d = density(sta[dnasig[ind]<8], adjust=.1)
    if(i==1) {
      plot(mylf(sta, at[ind]), type="l", main=paste("AT content (", chr, ")", sep=""),
           xlab="", ylab="", col="orange")
      plot(d, main=paste("Dim probes (", chr, ")", sep=""))
    } else {
      lines(d, col="blue")
    }
    
  }
}
cat("\n")
if(doPDF)
  dev.off()

## This code produces density of "eliminated" probes

if(FALSE) {
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
}
