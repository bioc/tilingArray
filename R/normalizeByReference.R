normalizeByReference = function(x, reference, whichBackground, nrStrata=20,
  cutoffQuantile=0.05, plotFileNames) {

  if(!is(x, "eSet"))
    stop("'x' must be an object of class 'eSet'")
  if(!is(reference, "eSet"))
    stop("'reference' must be an object of class 'eSet'")
  n = nrow(exprs(x))
  d = ncol(exprs(x))
  if(n!=nrow(exprs(reference)))
    stop("'x' and 'reference' must have 'exprs' elements in their 'assayData' slot with the same number of rows.")
  if(d<1)
    stop("There is nothing to normalize in 'x'.")
    
  ## reference normalization:
  refSig    = rowMeans(log(exprs(reference), 2))

  ## quantiles of the reference intensities, to group probes into
  ## strata for the background estimations
  quants    = quantile(refSig, probs=seq(0, 1, length=nrStrata+1))
  quants[1] = quants[1]-1

  ## reference signal just for the "background" probes
  refSigWhBg = refSig[whichBackground]

  ## strata is now a factor with 'nrStrata' levels and of same length as 'whichBackground'
  strata     = cut(refSigWhBg, quants)

  if(any(table(strata) < 5e3))
    warning("'some strata of background probes contain fewer than 5000 probes, are you sure this is allright?")
  
  xbg = (quants[-1]+quants[-length(quants)])/2  ## midpoint between quantiles
  ybg = matrix(as.numeric(NA), nrow=nrStrata, ncol=d)
  bgfun  = vector(mode="list", length=d)

  ## interpolate  
  for(j in 1:d) {
    ybg[, j] = tapply(log(exprs(x)[whichBackground, j], 2), strata, shorth)
    ## bgfun[[j]]  = locfit.raw(xbg, ybg[,j])
    bgfun[[j]] = approxfun(xbg, ybg[,j], rule=2)
  }

  ## diagnostic plot (also for the paper)
  if(!missing(plotFileNames)) {
    if(length(plotFileNames)!=d)
      stop("Please supply as many elements of 'plotFileNames' as there are arrays in 'x'")
    rgx = range(refSigWhBg)
    px  = seq(rgx[1], rgx[2], length=120)
    for(j in 1:d) {
      pdf(file=plotFileNames[j], width=8, height=6)
      smoothScatter(refSigWhBg, log(exprs(x)[whichBackground, j],2),
            xlab = "Reference intensity",
            ylab = "Background intensity", nrpoints=0)
      ## lines(px, predict(bgfun[[j]], newdata=px), col="darkred")
      lines(px, bgfun[[j]](px), col="darkred")
      dev.off()
    }
  }

  
  ## apply the background and the scaling
  cat("Applying background and scaling\n")
  yn = matrix(NA, nrow=nrow(exprs(x)), ncol=d)
  ttrefsig = 2^refSig
  for(j in 1:d) {
    ## yn[, j] = (exprs(x)[, j] - 2^predict(bgfun[[j]], newdata=refSig)) / ttrefsig
    yn[, j] = (exprs(x)[, j] - 2^bgfun[[j]](refSig)) / ttrefsig
  }
  
  e = new.env()
  if(d>=2) {
    cat("Between array normalization and variance stabilizing transformation\n")
    save(yn, file="yn.rda")
    browser()
    vsnres = vsn(yn, lts.quantile=0.95, subsample=2e5, niter=3) ## , verbose=FALSE
    assign("exprs", exprs(vsnres)/log(2), e)
    rm(vsnres)
  } else {
    warning("'x' has only one column, cannot do between array normalization and variance stabilizing transformation")
    assign("exprs", yn, e)
  }
  
  new("eSet", assayData=e, phenoData=phenoData(x), sampleNames=sampleNames(x))
}
