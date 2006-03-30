plotPenLL = function(seg, extrabar=numeric(0)) {
  cols = brewer.pal(3, "Set1")
  lls = cbind(logLik(seg), logLik(seg, penalty="AIC"), logLik(seg, penalty="BIC"))

  ylim = c(min(apply(lls, 2, quantile, probs=0.2)), max(lls))
  matplot(lls, type="b", col=cols, ylab="(penalized) log likelihood", lty=1, pch=16,xlab="S", ylim=ylim)

  optBIC = mean(which.max(lls[,3]))
  
  x0 = c(optBIC, extrabar)
  names(x0)[1] = cols[3]
  segments(x0=x0, x1=x0, y0=ylim[1]+diff(ylim)*0.05, y1=ylim[1]+diff(ylim)*0.95,
          col=names(x0), lty=2, lwd=2)
  legend(0, ylim[2], expression(log~L, log~tilde(L)[AIC], log~tilde(L)[BIC]),
         xjust=0, yjust=1, col=cols, lwd=1, pch=16)
}
