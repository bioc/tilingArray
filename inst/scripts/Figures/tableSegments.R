library("tilingArray")
source("readSegments.R") 

categorizeSegments = function(name) {
  stopifnot(!any(is.na(name)))
  
  care = c("gene (verif.)", "gene (unchar.)", "gene (dubious)",
           "other ncRNA", "unannotated")

  mt = matrix(0, nrow=length(name), ncol=length(care))
  colnames(mt)=care

  ## convert into regular expression
  pat = strsplit(name,  split=", ")
  pat[listLen(pat)==0] = "Unannotated"
  
  for(f in care) {
    theseNames = switch(f,
      "other ncRNA"    = gff$Name[gff$feature %in% c("ncRNA","snoRNA","snRNA", "tRNA", "rRNA")],
      ## orf_classification: applies to "gene", "CDS", and "intron"
      "gene (verif.)" = gff$Name[gff$feature=="gene" & gff$orf_classification=="Verified"],
      "gene (unchar.)" = gff$Name[gff$feature=="gene" & gff$orf_classification=="Uncharacterized"],
      "gene (dubious)" = gff$Name[gff$feature=="gene" & gff$orf_classification=="Dubious"],
      "unannotated"    = "Unannotated",
      stop("Zapperlot")
      )
    stopifnot(!is.null(theseNames))

    mt[, f] = sapply(pat, function(x) 
      !all(is.na(match(x, theseNames))))
  }
  ## clean up mt - set to zero all lower priority entries:
  for(i in 1:(ncol(mt)-1))
    mt[ mt[,i]>0, (i+1):ncol(mt) ] = 0 
    
  stopifnot(all(rowSums(mt)==1))
  mt
}


if(!exists("mt")) {
  graphics.off(); x11(width=9, height=5)
  par(mfrow=c(2,1))
  
  countDisjoint =function(sel) {sum(diff(sel)>1)}
  minOverlap=0.8
  maxDuplicated=0.5
  
  mt = vector(mode="list", length=length(rnaTypes))
  names(mt)=rnaTypes
  
  for(rt in rnaTypes) {
    s = get("segScore", get(rt))
    select = (s$frac.dup < maxDuplicated)
    unanno = (s$same.feature=="" | s$same.overlap < minOverlap)  & select
    thresh = calcThreshold(s$level, sel=unanno, showPlot=TRUE, main=rt)
    
    cat("=========", rt, "==========\n")
    cat("thresh=", signif(thresh, 2), "\n")

    transcribed = (select & (s$level>=thresh))
    mt[[rt]] = categorizeSegments(s$same.feature[transcribed])
  }

  dev.copy(pdf, "tableSegments-thresh.pdf", width=14, height=6); dev.off()
}

par(mfrow=c(1,2))
out = file("tableSegments.txt", open="wt")
for (i in seq(along=mt)) {
  cat("=========", rnaTypes[i], "==========\n", file=out)
  stopifnot(all(mt[[i]] %in% 0:1))
  freq = colSums(mt[[i]])
  nr   = nrow(mt[[i]])
  for(f in colnames(mt[[i]])) {
    cat(sprintf("%24s: %5d (%4.1f)\n", f,
                as.integer(freq[f]),
                100*freq[f]/nr), file=out)
  }
  cat(sprintf("%24s: %5d (%4.1f)\n", "total",
                nr, 100), file=out)

  ord = c(4, 1:3, 5)
  pie(freq[ord], main=c(polyA="poly-A RNA", tot="total RNA")[names(mt)[i]],
      col=c(brewer.pal(9, "Pastel1")[c(2:5, 1)])[ord],
      labels=paste(names(freq), " (", freq, ")", sep="")[ord])
}
close(out)
dev.copy(pdf, "tableSegments-pie.pdf", width=10, height=6); dev.off()

## x11()
## px = as.numeric(names(ss))
## plot(px, nrNew, ylim=c(0, max(nrNew, nrNewDj)), xlim=c(0, max(px)), type="b", pch=16)
## lines(px, nrNewDj, type="b", pch=16, col="blue")


