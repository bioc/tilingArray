library(tilingArray)

if(!exists("segwi"))load("1.+.rda")

## AIC = -2 logLik + 2*npar
## BIC = -2 logLik + log(ndat)*npar

penLL = function(ll, npar, ndat, what)
  switch(what,
  AIC = ll - npar,
  BIC = ll - log(ndat)*npar/2,
  stop("Zapperlot"))

plotPenLL = function(s) {
  cols = brewer.pal(3, "Set1")
  ll = logLik(s)
  plot(ll, type="l", col=cols[1], ylab="(penalized) log likelihood", xlab="S")
  lines(penLL(ll, npar=2*(1:length(ll)), ndat=nrow(s@y), what="AIC"), col=cols[2])
  lines(penLL(ll, npar=2*(1:length(ll)), ndat=nrow(s@y), what="BIC"), col=cols[3])
}


stop()

## using one segment: concordance of logLik.lm, my pedestrian formula,
## and the output of 'segment'

n  = seq(10, 100, by=10)
ll = matrix(as.numeric(NA), nrow=length(n), ncol=3)

myLogLik = function(res, n) {
  -n/2*( log(2*pi)+1-log(n)+log(sum(res*res)) )
}

for(i in seq(along=n)) {
  y = rnorm(n[i])
  lmy = lm(y~1)
  seg = segment(y, maxseg=2, maxk=length(y))
  res = y - coef(lmy)[1]
  ll[i, ] = c(stats::logLik(lmy), logLik(seg)[1], myLogLik(res, length(y)))
}

pairs(ll)
stopifnot(all(diff(t(ll))<1e-8))

## using multiple segments
 seg = segment(y, maxseg=10, maxk=length(y))
