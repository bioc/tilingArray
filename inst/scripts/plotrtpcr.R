library("RColorBrewer")
library(colorspace)
source("setScriptsDir.R")

graphics.off()
options(error=recover, warn=2)

rnaTypes = c("seg-tot-050909", "seg-polyA-050909", "seg-odT-050909")[3:2]
# 1. seg-tot-x is total-cDNA, while Lior's total-RH6 seems to refer to the total-RNA-RH6 samples;
# 2. Lior's poly-A samples seem to be the our so-called oligo-dT samples

source(scriptsDir("readSegments.R"))
source(scriptsDir("calcThreshold.R"))

gffOnlyGenes <- gff[gff$feature=="gene",]

if(!exists("rtpcrRes")) {
  rtpcrRes = read.table("Archi_RT_positions.txt", header=TRUE, as.is=TRUE, sep="\t")
  sp = strsplit(rtpcrRes$Locus, "_")
  gene = factor(sapply(sp, "[", 1))
  geneStrand <- sapply(levels(gene), function(x){
    x <- paste("^",x,"$",sep="")
    gffIdx <- union(grep(x, gffOnlyGenes$Name), grep(x, gffOnlyGenes$gene))
    if (length(gffIdx) != 1) warning(paste("Could not find unique Position for gene",x,"!"))
    return(as.character(gffOnlyGenes[gffIdx,"strand"]))})
  
  #rtpcrRes$arraytot = rtpcrRes$arrayodT = numeric(nrow(rtpcrRes))
  rtpcrRes$array1 <- numeric(nrow(rtpcrRes))

  for(j in 1:nrow(rtpcrRes)) {
    chr = rtpcrRes$Chr.[j]
    sta = rtpcrRes$strat[j]  ## sic
    end = rtpcrRes$end[j]
    
    myFun = function(dat) {
      sel = (dat$start>=sta) & (dat$end<=end) & (dat$unique==0)
      mean(dat$y[sel, ]) # note that this averages over replicates and probes
    }

    rtpcrRes$array1[j] = - get("theThreshold", get(rnaTypes[1])) + #0.5(
      myFun(get(paste(chr, geneStrand[gene[j]], "dat", sep="."), get(rnaTypes[1])))# +
      #myFun(get(paste(chr, "-", "dat", sep="."), get(rnaTypes[1])))
      # no averaging over strands (?): 
  }#for(j in 1:nrow(rtpcrRes))
}

## colors = rainbow(10)
colors = hex(HSV(seq(0, 360, length=11)[-11], 1, 1))

myPlot = function(x, y, fac, txtLabel,mypch=16, ...){
  stopifnot(length(x)==length(y), length(x)==length(fac))
  stopifnot(nlevels(fac)==length(colors))
  mycolors <- colors[as.integer(fac)]
  plot(x, y, col=mycolors, pch=mypch, ...)
  #plot(x, y, type="n",  ...)
  text(x, y, txtLabel, adj=c(0.5,0.5),col=mycolors, cex=1, pos=2)
  sp = split(seq(along=fac), fac)
  #for(j in seq(along=sp))
  #  lines(x[sp[[j]]], y[sp[[j]]], col=colors[j], lwd=2)
  ##browser()
  legend(min(x), max(y), levels(fac), col=colors, pch=16, cex=1, y.intersp=1)
  xyCor <- round(cor(x, y, method="spearman"),digits=2)
  text(max(x),min(y),labels=paste("Spearman Correlation =",xyCor), pos=2)
  # pos=2 means: put text to the left of the specified coordinates
}

graphics.off()
x11(width=10, height=10)
par(mfrow=c(1,1))

liorAnno = sapply(sp, "[", 2)
# according to Lior's one-to-two-letter annotation the first letter denotes
#  the array expression level of the PCR-amplified region, while the optional
#  second letter denotes the part of the gene in the case of segmented CDSs,
#  the so-called novel architectures.
exprLevel = substr(liorAnno, start=1, stop=1) 
region = substr(liorAnno, start=2, stop=2)

# plot symbol according to expr level defined by Lior:
mypch <- gsub("L", "25", as.character(exprLevel)) # triangle down
mypch <- gsub("H", "24", mypch) # triangle up
mypch <- as.numeric(gsub("M", "20", mypch))

myPlot(rtpcrRes$array1, -rtpcrRes$Mean.dT, gene, region,mypch=16, xlab="Array Expression (polyA)", ylab="rtPCR (polyA) [-Ct]")
# inverse ("-") of Mean.dT, since the number is, most likely, the PCR's Ct number, which denotes
#  the cycle number, when the amplified DNA reaches a certain amount. A higher DNA concentration in
#  the beginning thus relates to a lower Ct value (fewer amplifications necessary to reach detection
#  threshold). Therefore, the inverse should, ideally, yield a positive correlation in the plot.

#dev.copy(pdf, file="fig_rtcpr_b.pdf", width=10, height=10)
#dev.off()
