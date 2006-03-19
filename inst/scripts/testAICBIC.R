library(tilingArray)

if(!exists("segBig"))load("~/segBig.RData")

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
  ylim = ll[c(50, length(ll))]
  plot(ll, type="l", col=cols[1], ylab="(penalized) log likelihood", xlab="S", ylim=ylim)
  lines(penLL(ll, npar=2*(1:length(ll)), ndat=nrow(s@y), what="AIC"), col=cols[2])
  bic = penLL(ll, npar=2*(1:length(ll)), ndat=nrow(s@y), what="BIC")
  lines(bic, col=cols[3])

  lenChr= s@x[length(s@x)]
  nrSegChoice = round(lenChr/1500)
  optBIC = mean(which.max(bic))
  x0= c(nrSegChoice, optBIC)
  segments(x0=x0, x1=x0, y0=ylim[1]+diff(ylim)*0.05, y1=ylim[1]+diff(ylim)*0.95,
          col=c("#707070", cols[3]), lty=2)
  
  legend(0, ylim[2], expression(LL, tilde(LL)[AIC], tilde(LL)[BIC]),
         xjust=0, yjust=1, col=cols, lwd=1)
}

pdf(file="testAICBIC.pdf", width=4.5, height=4.5)
par(mai=c(1,1,0.01,0.01))
plotPenLL(segBig)
dev.off()

stop()

## Test the concordance of logLik.lm, my pedestrian formula,
## and the output of 'segment', for the trivial case of only one segment

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

