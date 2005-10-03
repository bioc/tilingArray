calcThreshold = function(x, sel, FDRthresh, dir, main) {
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

  if(!missing(dir)) {
    pdf(file=file.path(dir, "threshold2.pdf"), width=8, height=6)
  
    adjust = 0.5
    d1 = density(x[!sel], na.rm=TRUE, n=128, adjust=adjust)
    n  = length(d1$x)
    d2 = density(levu, na.rm=TRUE, from=d1$x[1], to=d1$x[n], n=n, , adjust=adjust)
    plot(d2, col="grey", lwd=2, main=main)
    lines(d1, col="red", lwd=2)
    abline(v=c(loc, thresh), col=c("black", "blue"))
    dn = dnorm(x=d1$x, mean=loc, sd=scale)
    lines(d1$x, dn/max(dn)*max(d2$y), col="orange")

    dev.off()
  }
  
  cat(main, ": loc=", signif(loc,4), "scale=", signif(scale,4), "thresh=", signif(thresh,4), "\n")
  return(thresh)
}



cat("Calculation of thresholds:\n",
    "==========================\n", sep="")

maxDuplicated = 0.5
FDRthresh     = 1e-3
cat("FDRthresh=", FDRthresh, "\n")

for(rt in rnaTypes) {

  alreadyDone = ("theThreshold" %in% ls(get(rt)))
  if(alreadyDone) {
    cat(rt, ": skipping threshold calculation since it was already done.\n", sep="")
  } else {
    s   = get("segScore", get(rt))
    sel = (s[, "frac.dup"] < maxDuplicated) & (s[, "overlappingFeature"] == "")
    thr = calcThreshold(s[, "level"], sel = sel, main=rt, FDRthresh=FDRthresh, dir=rt)
    s$level              = s[, "level"] - thr
    s$oppositeExpression = s[, "oppositeExpression"] - thr
    assign("segScore",  s,   envir=get(rt))
    assign("theThreshold", thr, envir=get(rt))
  }
}
