library("tilingArray")
source("scripts/readSegments.R") 
source("scripts/calcThreshold.R") 
source("scripts/categorizeSegments.R") 

options(error=recover, warn=0)
## out = file("tableSegments.txt", open="wt")
out = stdout()

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


colors = c(brewer.pal(9, "Pastel1")[c(2:6, 1)])
names(colors) =c("verified", "ncRNA", "uncharacterized", "dubious",
     "unA", "unI")

##
## PIE
##
if(0){
  par(mfrow=c(1,2))
  sink(out)
  cat("Counting the found segments:\n",
      "============================\n", sep="")
  for(rt in rnaTypes) {
    ct  = tab[[rt]]$count
    
    perct = 100*ct[,1]/ct[,2]
    ctp  = cbind(ct, "percent" = round(perct,1))
    ctp  = rbind(ctp,  "total" = colSums(ctp, na.rm=TRUE))
    
    cat("\n", rt, "\n-----\n", sep="")
    print(ctp)
    
    px = ct[, "observed"]
    names(px) = rownames(ct)
    px = px[c(1,3,2,4,6,5)]
    pie(px, radius=0.75, main=c(polyA="poly-A RNA", tot="total RNA")[rt],
        col=colors,
        labels=paste(names(px), " (", px, ")", sep=""))
  }
  sink()
  dev.copy(pdf, "tableSegments-pie.pdf", width=14, height=4.8); dev.off()
}

##
## LENGTH DISTRIBUTIONS
##
if(0){
par(mfrow=c(2,5))
maxlen=6000
br = seq(0, maxlen, by=100)
for(rt in rnaTypes) {
  s  = get("segScore", get(rt))
  ct = tab[[rt]]$category
  for(lev in levels(ct)[c(1:3, 5:6)]) {
    len = s$length[ct == lev]
    len[len>maxlen]=maxlen
    hist(len, breaks=br, col=colors[lev], main=paste(rt, lev))
  }
  
}
dev.copy(pdf, "tableSegments-lengths.pdf", width=14, height=6); dev.off()
}

##
## CONSERVATION
##

blastResultFiles = c("Sbay_contigs.out", "Smik_contigs.out", "Spar_contigs.out", "Spom_all.out")
names(blastResultFiles) = c("S.bayanus", "S.mikatae", "S.paradoxus", "S.pombe")
if(!exists("blastres")) {
  blastres = lapply(blastResultFiles, function(f)
	 read.table(file.path(indir, "fasta", f),
              sep="\t", as.is=TRUE, header=FALSE))
}	
sink(out)
cat("\n\nPercent of segments with BLAST-hit in other genome (E<",
    ethresh, "):\n", "====================================================\n", sep="")


for(rt in rnaTypes[1]) {
  s  = get("segScore", get(rt))
  ct = tab[[rt]]$category
  
  sp = split(seq(along=ct), ct)

  levthresh = median(s$level, na.rm=TRUE)
  sp = lapply(sp, function(i) i[s$level[i]>levthresh])
  
  sp[[7]] = seq(along=ct)
  names(sp)[7] = "whole genome"
  
  hit = matrix(NA, nrow=length(sp), ncol=length(blastResultFiles))
  ethresh = 1e-20
  rownames(hit) = names(sp)
  colnames(hit) = names(blastResultFiles)
  for(b in seq(along=blastres)) {
    br = blastres[[b]]
    br = br[br[[11]]<ethresh, ]
    hit[,b] = sapply(sp, function(x) mean(x %in% br[[1]]))
  }
  cat("\n", rt, "\n-----\n", sep="")
  print(signif(hit,2))
}
  sink()


is(!is(out, "terminal"))
  close(out)


