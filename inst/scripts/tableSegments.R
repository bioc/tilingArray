library("tilingArray")
source("scripts/readSegments.R") 
source("scripts/calcThreshold.R") 
source("scripts/categorizeSegments.R") 

if(!exists("tab")) {
  graphics.off(); x11(width=9, height=5)
  par(mfrow=c(2,1))
  
  tab = vector(mode="list", length=length(rnaTypes))
  names(tab)=rnaTypes
  
  for(rt in rnaTypes) {
    s = get("segScore", get(rt))
    tab[[rt]] = categorizeSegmentsPie(s)
  }

  dev.copy(pdf, "tableSegments-thresh.pdf", width=11, height=8); dev.off()
}

par(mfrow=c(1,2))
out = file("tableSegments.txt", open="wt")
for (i in seq(along=tab)) {
  cat("=========", rnaTypes[i], "==========\n", file=out)
  mt = tab[[i]]$mt
  
  stopifnot(all(mt %in% 0:1))
  freq = colSums(mt)
  nr   = nrow(mt)
  for(f in colnames(mt)) {
    cat(sprintf("%14s: %5d (%4.1f percent) -- in genome: %5d\n", f,
                as.integer(freq[f]),
                100*freq[f]/nr, tab[[i]]$nrInGenome[f]), file=out)
  }
  cat(sprintf("%14s: %5d (%4.1f)\n", "total",
                nr, 100), file=out)

  ord = c(4, 1:3, 5)
  pie(freq[ord], radius=0.75, main=c(polyA="poly-A RNA", tot="total RNA")[names(tab)[i]],
      col=c(brewer.pal(9, "Pastel1")[c(2:5, 1)])[ord],
      labels=paste(names(freq), " (", freq, ")", sep="")[ord])
}
close(out)
dev.copy(pdf, "tableSegments-pie.pdf", width=14, height=4.8); dev.off()

## x11()
## px = as.numeric(names(ss))
## plot(px, nrNew, ylim=c(0, max(nrNew, nrNewDj)), xlim=c(0, max(px)), type="b", pch=16)
## lines(px, nrNewDj, type="b", pch=16, col="blue")


