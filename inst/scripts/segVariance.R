options(error=recover, warn=2)
library("tilingArray")
library("prada")

indir = "segmentation-3polyA"
nrchr = 3
js = paste(rep(1:nrchr, each=2), rep(c("+", "-"), nrchr), "rda", sep=".")

if(!exists("s")) {
  s  = new.env()
  r  = new.env()
  cat("Loading.\n")
  for(j in js) {
    load(file.path(indir, j))
    assign(sub("rda", "seg", j), seg, envir=s)
    assign(sub("rda", "dat", j), dat, envir=s)
    load(file.path(indir, sub(".rda", "rand.rda", j)))
    assign(sub("rda", "seg", j), seg, envir=r)
    assign(sub("rda", "dat", j), dat, envir=r)
  }
}

calcStat = function(e, nrBasePerSeg = 1500) {
  dat = get(sub("rda", "dat", j), envir=e)
  seg = get(sub("rda", "seg", j), envir=e)
  stopifnot(nrow(dat$y)==length(dat$x))
  cp = round(max(dat$x)/nrBasePerSeg)
  th = c(1, seg$th[cp, 1:cp])
  ct = cut(1:nrow(dat$y), breaks=th-1)
  ## note: dat$y is a matrix with length(ct) rows. The recycling rule
  ## works in our favour here!
  sp = split(dat$y, ct)
  v  = sapply(sp, var)
  m  = sapply(sp, mean)
  n  = listLen(sp)
  
  if(TRUE){
    t = numeric(length(sp))
    for (z in 2:length(sp)) {
      sdev = sqrt(  ((n[z-1]-1)*v[z-1]+(n[z]-1)*v[z]) / (n[z]+n[z-1]-2) )
      t[z] = (m[z]-m[z-1]) / (sdev*sqrt(1/n[z]+1/n[z-1]))
    }
    z = sample(2:length(sp), 1)
    stopifnot(abs(t[z] - t.test(sp[[z]], sp[[z-1]], var.equal=TRUE)$statistic) < 1e-9)
  }

  if(FALSE)
    t=c(0, diff(m))
  
  return(list(statistic=t, n=n))
}

graphics.off(); x11()
j = js[3]

segs = get(sub("rda", "seg", j), envir=s)
segr = get(sub("rda", "seg", j), envir=r)
par(mfrow=c(2,2))
plot(segr$J, segs$J)
plot(2*seq(along=segs$J), segs$J, main="data")
plot(2*seq(along=segr$J), segr$J, main="random")
ddJ = diff(segs$J) ## diff(diff(segs$J))
plot(2*seq(along=ddJ), ddJ, ylim=quantile(ddJ, c(0.01, .99), na.rm=TRUE), main="data")
stop()


ts = calcStat(s)
tr = calcStat(r)


par(mfrow=c(3,2))
hist(ts$statistic, col="lightblue", breaks=50)
hist(tr$statistic, col="orange", breaks=50)
hist(ts$n, col="lightblue", breaks=50)
hist(tr$n, col="orange", breaks=50)
qqplot(ts$statistic, tr$statistic); abline(a=0, b=1, col="red")


if(FALSE) {
  v = sapply(js, function(j)
     sd(get(sub("rda", "dat", j), envir=s)$y))
}

## see whether variance depends on position along chromosome
## (it doesn't seem to)
if(FALSE) {
  par(mfrow=c(4,4), mai=c(0.3, 0.4, 0.4, 0.01))
  for(j in js[17:32]) {
    dat = get(sub("rda", "dat", j), envir=s)
    dy  = diff(dat$y)
    dy  = dy[abs(dy)<=2]
    smoothScatter(seq(along=dy), dy, pch=".", main=j)
  }
}
