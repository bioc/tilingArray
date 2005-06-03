

calcThreshold = function(x, sel, FDRthresh, showPlot=FALSE, main) {
  stopifnot(is.numeric(x), is.logical(sel))

  require("multtest")
  require("genefilter")

  levu  = x[sel]
  loc   = shorth(levu, na.rm=TRUE)
  z     = levu[which(levu<=loc)]-loc
  scale = mad(c(z, -z))
  ## we calculate the threshold for level on the basis of the normal
  ## distribution and the FDR for _unannotated_ features
  pn = pnorm(q=x, mean=loc, sd=scale, lower.tail=FALSE)
  bh = mt.rawp2adjp(pn[sel], proc="BY")
  stopifnot(all(bh$adjp[, 1] == pn[sel][bh$index], na.rm=TRUE))
  
  adjp = numeric(nrow(bh$adjp))
  adjp[bh$index] = bh$adjp[,2]
  selfdr = (adjp < FDRthresh)

  thresh = min(x[sel][selfdr], na.rm=TRUE)

  if(showPlot) {
    adjust = 0.5
    d1 = density(x[!sel], na.rm=TRUE, n=128, adjust=adjust)
    n  = length(d1$x)
    d2 = density(levu, na.rm=TRUE, from=d1$x[1], to=d1$x[n], n=n, , adjust=adjust)
    plot(d2, col="grey", lwd=2, main=main)
    lines(d1, col="red", lwd=2)
    abline(v=c(loc, thresh), col=c("black", "blue"))
    dn = dnorm(x=d1$x, mean=loc, sd=scale)
    lines(d1$x, dn/max(dn)*max(d2$y), col="orange")
  }
  cat(main, ": loc=", signif(loc,3), "scale=", signif(scale,3), "thresh=", signif(thresh,3), "\n")
  browser()
  return(thresh)
}

cat("Calculating Thresholds:\n",
    "=======================\n", sep="")

if(interact) {
  ## x11(width=10, height=length(rnaTypes)*3)
} else {
  pdf(file="tableSegments-thresh.pdf", width=11, height=length(rnaTypes)*4)
}
par(mfrow=c(length(rnaTypes),1))

maxDuplicated = 0.5
FDRthresh     = 1e-3
cat("FDRthresh=", FDRthresh, "\n")

for(rt in rnaTypes) {
  s = get("segScore", get(rt))
  sel = (s[, "frac.dup"] < maxDuplicated) & (s[, "overlappingFeature"] == "")
  thr = calcThreshold(s[, "level"], sel = sel, main=rt, FDRthresh=FDRthresh, showPlot=TRUE) ##  showPlot=!interact)
  assign("threshold", thr, envir=get(rt))
}

if(!interact)
    dev.off()

cat("\n\n")
