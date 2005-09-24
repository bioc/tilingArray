##
## run this script after "tableSegments"
##
library("GOstats")
source(scriptsDir("GOHyperG.R"))

interact= !TRUE

outfile = "antisense-GO"
if(!interact)
  sink(paste(outfile, "txt", sep="."))

asCand = matrix(FALSE, nrow=length(featNames$"annotated ORFs"), ncol=2)
rownames(asCand)=featNames$"annotated ORFs"
colnames(asCand)=c("filtered", "all")

for(what in colnames(asCand)) {
  catgSel = list(filtered = "novel antisense - filtered",
           all = c("novel antisense - filtered", "novel antisense - unassigned"))[[what]]

  count = numeric(2)
  hit   = numeric(3)
  names(hit) = c("5'", "3'", "tot")
  
  cat("Analysis of level differences:\n")
  cat("=====================================\n")
  for(rt in rnaTypes) {
    s = cs[[rt]]
    selSeg = which(s[, "category"] %in% catgSel)
    asGenes = unique(unlist(strsplit(s[selSeg, "oppositeFeature"], split=", ")))
    asGenes = intersect(asGenes, rownames(asCand))
    asCand[ asGenes, what ] = TRUE

    ## count how often antisense overlaps 3' and 5' ends of gene
    isGeneSegment = (s$featureInSegment %in% featNames$"annotated ORFs") & (!is.na(s$level))
    isAntiSenseGeneSegment = rep(FALSE, length(isGeneSegment))
  
    for(j in selSeg) {
      of = strsplit(s[j, "oppositeFeature"], ", ")[[1]]
      stopifnot(length(of)>=1)
      if(length(of)==1) {
        w = which(gff$Name==of & isGene)
        stopifnot(length(w)<=1)
        if(length(w)==1) {
          delta   = c((s[j, "start"]-50)<=gff$start[w],
                      (s[j, "end"]+50)  >=gff$end[w],
            TRUE) 
          strand = as.character(gff$strand[w])
          stopifnot(strand%in%c("+", "-"))
          if(strand=="-") ## reverse
            delta[1:2]=delta[2:1]
          hit = hit + delta
          
          v = which(s$featureInSegment==of)
          stopifnot(all(isGeneSegment[v]))
          isAntiSenseGeneSegment[v] = TRUE
          
        } else {
          count[2]=count[2]+1
        }
      } else {
        count[1]=count[1]+1
      }
    }
    x1 = s$level[isGeneSegment & !isAntiSenseGeneSegment]
    x2 = s$level[isGeneSegment &  isAntiSenseGeneSegment]
    cat(rt, " & ", catgSel, ": ",
        "median level: ", signif(median(x2), 3), " vs ",
        signif(median(x1), 3), "   p=", format.pval(wilcox.test(x1, x2)$p.value), "\n\n", sep="")
  } ## for rt
  
  cat(sum(asCand[, what]), " genes were found opposite a segment of category ",
      paste("'", catgSel, "'", sep="", collapse=", "), ".\n\n", sep="")

  cat("Analysis of 3'/5' preference of potential antisense transcripts:\n")
  cat("================================================================\n")
  cat("Number of cases where more than one feature in oppositeFeature:", count[1], "\n")
  cat("Number of cases where oppositeFeature not a gene:", count[2], "\n")

  cat("hit:\n")
  print(hit)
  cat("\n\n\n")

} ## for what

##--------------------------------------------------
## Do antisense genes have longer/shorter UTRs?
##--------------------------------------------------
cat("\n\nDo antisense genes have longer/shorter UTRs?\n")
cat("================================================\n\n")
if(!exists("utr"))
  load("utr.rda")
utrls = utr[["combined"]]
for(j in colnames(utrls)) {
  for(what in colnames(asCand)) {
    x1 = utrls[ (rownames(utrls) %in% names(which(asCand[, what]))), j] 
    x2 = utrls[!(rownames(utrls) %in% names(which(asCand[, what]))), j] 
    p = wilcox.test(x1, x2)$p.value
    cat(j, " & ", what, ": median length=", median(x1), " vs ", median(x2), "   p=", format.pval(p), "\n", sep="")
  }
}


what="all"
cat("\n\nGO category analysis (", what, "):\n", sep="")
cat("======================================\n\n")
GOHyperG(rownames(asCand)[asCand[, what]])

if(!interact) {
  dev.copy(pdf, paste(outfile, what, "pdf", sep="."), width=12, height=6.3)
  dev.off()
}



if(!interact)
  sink()

