segChrom = function(y, probeAnno, chr=1:17, strands=c("+", "-"),
  nrBasesPerSegment = 1500, maxk = 3000, step = 7, confint = FALSE,
  confintLevel = 0.95,   useLocks=TRUE, verbose=TRUE, savedir) {
  
  segObj = new.env(parent = baseenv())
  chrstrd = paste(rep(chr, each=length(strands)), rep(strands, times=length(chr)), sep=".")
  
  for(j in seq(along=chrstrd)) {        

    skip=FALSE
    if(!missing(savedir)) {
      filename = file.path(savedir, paste(chrstrd[j], "rda", sep="."))
      if(useLocks)
        if(file.exists(filename)) {
          skip=TRUE
          if(verbose)
            cat(sprintf("Skipping %s since %s already exists.\n", chrstrd[j], filename))
        } else {
          con=file(filename, open="wt")
          cat(system("uname -a; echo \"\n\"; date", intern=TRUE), file=con)
          close(con)
        }
    }

    if(!skip){
      if(verbose)
        cat(sprintf("Running 'segment' on chromosome %s", chrstrd[j]))

      what = c("start", "end", "index", "unique")
      prbs = do.call("data.frame", mget(paste(chrstrd[j], what, sep = "."), probeAnno))
      colnames(prbs) = what
      prbs$mid = (prbs$start + prbs$end)/2
      prbs = prbs[order(prbs$mid), ]
      
      if(is.matrix(y))
        numna = rowSums(is.na(y[prbs$ind, ]))
      else
        numna = rowSums(is.na(exprs(y)[prbs$ind, ]))
      stopifnot(all(numna %in% c(0, ncol(y))))
    
      prbs = prbs[numna == 0, ]
      sprb = prbs[sampleStep(prbs$mid, step = 7), ]         
      nsegs = as.integer(round(sprb$end[nrow(sprb)]/nrBasesPerSegment))
    
      if(is.matrix(y))
        ychr = y[sprb$ind, ,drop = FALSE]
      else
        ychr = exprs(y)[sprb$ind, ,drop = FALSE]
    
      s = segment(ychr, maxseg = nsegs, maxk = maxk)
      
      s@x = sprb$mid
      s@flag = sprb$unique
      if(confint)
        assign(chrstrd[j], confint(s, parm = nsegs, level = confintLevel), segObj)
      else {
        s@nrSegments = nsegs
        assign(chrstrd[j], s,segObj)
      }
      
    if(!missing(savedir))
      save(s, file=filename)
      
    if(verbose)
      cat(" ... complete\n")
    }
  } ##  for j
  segObj
}
