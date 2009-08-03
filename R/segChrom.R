segChrom = function(y, probeAnno, chr=1:17, strands=c("+", "-"),
		nrBasesPerSegment = 1500, maxk = 3000, step = 7, confint = FALSE,
		confintLevel = 0.95,   useLocks=TRUE, verbose=TRUE, savedir) {
	
  segObj = new.env(parent = baseenv())
  chrstrd = if((length(strands)==1) && is.na(strands)) {
    chr 
  } else { 
    paste(rep(chr, each=length(strands)), rep(strands, times=length(chr)), sep=".")
  }
	
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
          ## create lock file
	  con=file(filename, open="wt")
	  cat(system("uname -a; echo \"\n\"; date", intern=TRUE), file=con)
	  close(con)
        }
    } ## if(useLocks)
    
    if(!skip){
      if(verbose)
        cat(sprintf("Running 'segment' on chromosome %s", chrstrd[j]))
      
      ## construct dataframe 'df' with four columns as follows		
      df_colnames = c("start", "end", "index", "unique")
      pa_elements = paste(chrstrd[j], df_colnames, sep = ".")
      prbs = if(is(probeAnno, "probeAnno")) {
        lapply(pa_elements, function(w) probeAnno[w])
      } else if (is(probeAnno, "environment")) {
        mget(pa_elements, probeAnno)
      } else {
        stop(sprintf("Invalid class of argument 'probeAnno', is '%s', should be class 'probeAnno' or 'environment'.", class(probeAnno)))  
      } 	     

      prbs = do.call(data.frame, prbs)
      colnames(prbs) = df_colnames
			
      ## sort probes by chromosomal midpoint position:
      prbs$mid = (prbs$start + prbs$end)/2
      prbs = prbs[order(prbs$mid), ]
      
      ## remove missing (NA) values:
      numna = rowSums(is.na((if(is.matrix(y)) y else exprs(y))[prbs$index,,drop = FALSE]))
      stopifnot(all(numna %in% c(0, ncol(y))))
      prbs = prbs[numna == 0, ]
			
      ## subsample probes to overcome irregular probe spacing
      ##  esp. dense spacing in repetitive regions:
      sprb = prbs[sampleStep(prbs$mid, step = step), ]
      
      ## determine number of segments
      nsegs = as.integer(round(sprb$end[nrow(sprb)]/nrBasesPerSegment))
			
      ychr = (if(is.matrix(y)) y else exprs(y))[sprb$index, ,drop = FALSE]
	
      s = segment(ychr, maxseg = nsegs, maxk = maxk)
			
      s@x = sprb$mid
      s@flag = as.integer(sprb$unique)
			
      if(confint) {
        assign(chrstrd[j], confint(s, parm = nsegs, level = confintLevel), segObj)
      } else {
        s@nrSegments = nsegs
        assign(chrstrd[j], s, segObj)
      }
			
      if(!missing(savedir))
        save(s, file=filename)
			
      if(verbose)
        cat(" ... complete\n")
    } ## if(!skip)
  } ##  for j	
  return(segObj)
}

