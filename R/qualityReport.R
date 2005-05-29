## Yeast Tiling Array Project, March 2005
## script to quickly check quality of provided new samples

## at present this script is meant to process new RNA samples from the
##  S96 strain.

## x: affyBatch or name of celfiles
qualityReport = function(x, hybeType, normRef=NULL, compress = TRUE,
                         selectGenes,
                         gff, probeAnno, 
                         output    = "HTML",
                         outputDir = ".",
                         verbose = TRUE)
{
  ## check parameter values and set up:
  switch(output,
    HTML =  {
      if(!require(R2HTML))
        stop("Could not load required package R2HTML!")
    },
    stop(paste("output format '", output, "' is not implemented.", sep=""))
  ) ## switch
  
  if(!file.exists(outputDir))
    stop(paste("outputDir '", outputDir, "' does not exist.", sep=""))
  if(!file.info(outputDir)$isdir)
    stop(paste("'", outputDir, "' is not a directory.", sep=""))

  graphicsDir = paste(outputDir,"png",sep="/")
  if(!file.exists(graphicsDir))
    dir.create(graphicsDir)
  
  if(is(x, "character"))
    x = read.affybatch(filenames=x, compress=compress, verbose=verbose)
  if(!is(x, "AffyBatch"))
    stop("'x' must be an AffyBatch")

  if(length(hybeType)!=1 || !(hybeType %in% c("Direct", "Reverse")))
    stop(paste("'hybeType' must have same length 1 and be",
               "either 'Direct' or 'Reverse'."))
  
  ## log2
  exprs(x) = log(exprs(x), 2)

  ## normalize if applicable
  if(hybeType=="Direct") {
    cat("--> Not normalizing since hybeType=Direct. <--\n")
    normRef = NULL
  }
  if(!is.null(normRef)) {
    if(is(normRef, "character"))
      normRef = read.affybatch(filenames=normRef, compress=compress, verbose=verbose)
    if(!is(normRef, "AffyBatch"))
      stop("'normRef' must be an AffyBatch")
    xu = x ## save unnormalized data
    exprs(x) = exprs(x) - log(rowMeans(exprs(normRef)), 2)
  }

  sampleNames = gsub("\\.cel(\\.gz)*$", "", rownames(pData(x)))
  sampleNames = gsub("#|&|/|\\*", "_", sampleNames)
  
  ## plot some values along the chromosome for selected ORFs:
  gff$Name = getAttributeField(gff$attributes, "Name")
  selG = which(gff$Name %in% selectGenes & gff$feature=="gene")
  stopifnot(length(selG)==length(selectGenes))
  
  ## work on individual samples:
  for (s in 1:nrow(pData(x))) {
    if (verbose)
      cat("Working on ", sampleNames[s], ":\n", sep="")
    if (output=="HTML") {
        tit = paste("Quality report for", sampleNames[s])
        out = HTMLInitFile(outdir = outputDir,
              filename = sampleNames[s],
              Title = tit, CSSFile="R2HTML.css")
        HTML.title(tit, HR=1, file=out)
        HTML.title("Along Chromosome Plots", HR=3, file=out)
      }

    par(mfrow=c(1, 1))
    for (i in seq(along=selG)) {
      wh = selG[i]
      plotAlongChrom(y    = exprs(x)[,s],
                     hybeType = hybeType,
                     chr  = gff$chr[wh], 
                     from = gff$start[wh],
                     to   = gff$end[wh],
                     extend=3000,
                     gff  = gff,
                     probeAnno = probeAnno)
      if(output == "HTML") {
        HTMLplot(file=out, Width=700, Height=300, GraphDirectory=graphicsDir,
                 GraphFileName=paste(sampleNames[s], selectGenes[i], sep="-"),
                 GraphBorder=0)
      }
    } ## for i

    ## calculate scores
    probe = get(paste("probe", hybeType, sep=""), probeAnno)
    nsc = calcScores(exprs(x), probe)

    ## compare with unnomalized score (if applicable)
    if(!is.null(normRef)) {
      usc = calcScores(exprs(xu), probe)

      myplot=function(x, y, ...) {
        axlim = quantile(c(x,y), c(0.01,0.99))
        plot(x, y, xlim=axlim, ylim=axlim, pch=".",
             xlab = "unnormalized", ylab = "normalized", ...)
        abline(a=0, b=1, col="red")
      }

      par(mfrow=c(1, 2))
      myplot(usc$sd, nsc$sd, main="standard deviation")
      p1 = -log(usc$p, 10)
      p2 = -log(nsc$p, 10)
      p1[p1>100]=100
      p2[p2>100]=100
      myplot(p1, p2, main="-log10(pvalue)")
      
      if (output == "HTML") {
        HTMLhr(file=out)
        HTML.title("Comparing CDSs between raw and normalized data", HR=3, file=out)
        HTMLplot(file=out, Width=800, Height=450, GraphDirectory=graphicsDir,
                 GraphFileName=paste(sampleNames[s], "norm", sep="-"),
                 GraphBorder=0)
      }
    }
    
    ## histogram of levels
    par(mfrow=c(1, 1))
    hist(nsc$mean, main="Histogram of mean levels", col="lightblue")

    pthresh = 0.05/length(nsc$p)
    numberExpr = paste(sum(nsc$p < pthresh), " CDSs have p < ", signif(pthresh,2),
        " =(0.05/", length(nsc$p), ")\n", sep="")

    if (output == "HTML"){
      HTMLhr(file=out)
      HTML.title("Histogram of CDSs' mean levels", HR=3, file=out)
      HTMLplot(file=out, Width=600, Height=400, GraphDirectory=graphicsDir,
               GraphFileName=paste(sampleNames[s], "hist", sep="-"),
               GraphBorder=0)
      HTML(numberExpr, file=out)
      HTMLEndFile(file=out)

    }
    
  } # for (s in 1:nsamples)
  invisible(NULL)
}


## t-scores and standard deviations for annotated CDSs
calcScores = function(x, probe) {
  browser()
  baseline = median(x[probe$no_feature=="no"])
  x = x-baseline
  xs    = split(x, probe$CDS)
  ## throw out CDS with less than 7 probes
  xs    = xs[names(xs)!="" & listLen(xs) > 7]
  lls   = listLen(xs)
  means = sapply(xs, mean)
  sds   = sapply(xs, sd)
  p     = pt(means/sds*sqrt(lls), df=lls-1, lower.tail=FALSE)
  list(mean = means, sd = sds, p = p)  
}
