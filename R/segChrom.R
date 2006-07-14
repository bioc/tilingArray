segChrom = function(y, probeAnno, chr=1, nsegs, nrBasesPerSegment = 1500, step=7, 
                    confint=FALSE, confintLevel=0.95, save=FALSE, savedir=NULL) {
     segObj = new.env(parent = baseenv())
     strand = c("+", "-")
     chrstrd = paste(rep(chr, each=2), rep(strand, times=length(chr)), sep=".")
     if(missing(nsegs))
       nsegs = rep(NA, length(chrstrd))
     for(j in seq(along=chrstrd)) {        
         cat(sprintf("Running segment on chromosome %s", chrstrd[j]))
         what = c("start", "end", "index", "unique")
         prbs = do.call("data.frame", mget(paste(chrstrd[j], what, sep = "."), probeAnno))
         colnames(prbs) = what
         prbs$mid = (prbs$start + prbs$end)/2
         prbs = prbs[order(prbs$mid), ]
         if(is.matrix(y))
           numna = rowSums(is.na(y[prbs$ind, ]))
         else
           numna = rowSums(is.na(exprs(y)[prbs$ind, ]))
#         stopifnot(all(numna %in% c(0, ncol(y))))
         prbs = prbs[numna == 0, ]
         sprb = prbs[sampleStep(prbs$mid, step = 7), ]         
         if(is.na(nsegs[j]))
           nsegs[j] = round(sprb$end[nrow(sprb)]/nrBasesPerSegment)
         if(is.matrix(y))
           ychr = y[sprb$ind, ,drop = FALSE]
         else
           ychr = exprs(y)[sprb$ind, ,drop = FALSE]
         s = segment(ychr, maxseg = nsegs[j], maxk = 3000)
         s@x = sprb$mid
         s@flag = sprb$unique
         if(confint)
           assign(chrstrd[j], confint(s, parm = as.integer(nsegs[j]), level = confintLevel), segObj)
         else {
           s@nrSegments = as.integer(nsegs[j])
           assign(chrstrd[j], s,segObj)
           }
         if(save)
           save(s, file=file.path(ifelse(is.null(savedir), ".", savedir), paste(chrstrd[j], "rda", sep=".")))
         cat(" ... complete\n")
         }
      segObj
      }
