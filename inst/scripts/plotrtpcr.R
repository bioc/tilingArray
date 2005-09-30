library("RColorBrewer")
source("setScriptsDir.R")

graphics.off()
options(error=recover, warn=2)

rnaTypes = c("seg-tot-050909", "seg-polyA-050909")[1]
## rnaTypes = c("seg-tot-050909", "seg-odT-050909")

source(scriptsDir("readSegments.R"))
source(scriptsDir("calcThreshold.R"))

if(!exists("rtpcrRes")) {
  rtpcrRes = read.table("Archi_RT_positions.txt", header=TRUE, as.is=TRUE, sep="\t")

  rtpcrRes$arraytot = rtpcrRes$arrayodT = numeric(nrow(rtpcrRes))

  for(j in 1:nrow(rtpcrRes)) {
    chr = rtpcrRes$Chr.[j]
    sta = rtpcrRes$strat[j]  ## sic
    end = rtpcrRes$end[j]
    
    myFun = function(dat) {
      sel = (dat$start>=sta) & (dat$end<=end) & (dat$unique==0)
      mean(dat$y[sel, ]) # note that this averages over replicates and probes
    }

    ## 
    rtpcrRes$array1[j] = - get("theThreshold", get(rnaTypes[1])) + 0.5*(
                       myFun(get(paste(chr, "+", "dat", sep="."), get(rnaTypes[1]))) +
                       myFun(get(paste(chr, "-", "dat", sep="."), get(rnaTypes[1]))))
    
  }
}

## colors = c(brewer.pal(11, "Set3")[-2])
colors = rainbow(10)

myPlot = function(x, y, fac, txtLabel, ...){
  stopifnot(length(x)==length(y), length(x)==length(fac))
  stopifnot(nlevels(fac)==length(colors))
  ## plot(x, y, col=colors[as.integer(fac)], pch=16, xlab=xlab, ylab=ylab, ...)
  plot(x, y, type="n",  ...)
  text(x, y, txtLabel, adj=c(0.5,0.5),col=colors[as.integer(fac)], cex=1)
  sp = split(seq(along=fac), fac)
  for(j in seq(along=sp))
    lines(x[sp[[j]]], y[sp[[j]]], col=colors[j], lwd=2)
  ##browser()
}

graphics.off()
x11(width=10, height=10)
par(mfrow=c(1,1))

sp = strsplit(rtpcrRes$Locus, "_")
gene = factor(sapply(sp, "[", 1))
region = sapply(sp, "[", 2)
  
myPlot(rtpcrRes$array1, rtpcrRes$Mean.RH6, gene, region,  xlab="array (total)", ylab="rtPCR (total+RH6)")
legend(-0.5, 19.5, levels(gene), col=colors, pch=16, cex=1, y.intersp=1)

dev.copy(pdf, file="fig_rtcpr_b.pdf", width=10, height=10)
dev.off()
