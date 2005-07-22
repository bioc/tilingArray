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

  graphicsDir = "png"
  if(!file.exists(file.path(outputDir, graphicsDir)))
    dir.create(file.path(outputDir, graphicsDir))
  
  if(is(x, "character"))
    x = read.affybatch(filenames=x, compress=compress, verbose=verbose)
  if(!is(x, "AffyBatch"))
    stop("'x' must be an AffyBatch")

  if(length(hybeType)!=1 || !(hybeType %in% c("Direct", "Reverse")))
    stop(paste("'hybeType' must have length 1 and be",
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
        HTMLplot(file=out, Width=700, Height=300, GraphDirectory=outputDir,
                 GraphFileName=file.path(graphicsDir, paste(sampleNames[s], selectGenes[i], sep="-")),
                 GraphBorder=0)
      }
    } ## for i

    ## calculate scores
    probe = get(paste("probe", hybeType, sep=""), probeAnno)
    ## nsc = calcScores(exprs(x), probe)


    if (output == "HTML"){
      HTMLEndFile(file=out)
    }
    
  } # for (s in 1:nsamples)
  invisible(NULL)
}
