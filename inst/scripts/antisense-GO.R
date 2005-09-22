##
## run this script after "tableSegments"
##
library("GOstats")
source(scriptsDir("GOHyperG.R"))

what = c("filtered", "all")[1]
switch(what,
       filtered = {
         outfile = "antisense-GO-filt"
         catgSel = "novel antisense - filtered"
       },
       all = {
         outfile = "antisense-GO-all"
         catgSel = c("novel antisense - filtered", "novel antisense - unassigned")
       },
       stop("Zapperlot"))


interact=!TRUE
if(!interact)
  sink(paste(outfile, "txt", sep="."))

asCand = NULL
count1 = count2 = count3 = 0
hit    = numeric(3)
names(hit) = c("5'", "3'", "tot")
isGene  = ((gff$feature=="gene") & (gff$orf_classification %in% c("Uncharacterized", "Verified")) &
     (gff$chr<=16))  
allGenes = unique(gff$Name[isGene])

  
for(rt in rnaTypes) {
  s = cs[[rt]]
  selSeg     = which(s[, "category"] %in% catgSel)
  asCand     = c(asCand, s[ selSeg, "oppositeFeature"])

  isGeneSegment = (s$featureInSegment %in% allGenes) & (!is.na(s$level))

  ## ovfs = strsplit(s$overlappingFeature, ", ")
  ## isGeneSegment = sapply(ovfs, function(x) any(x %in% allGenes)) & (!is.na(s$level))
  
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
        stopifnot(length(v)<=1,  isGeneSegment[v])
        if(length(v)==1) {
          isAntiSenseGeneSegment[v] = TRUE
        } else {
          count3=count3+1
        }      
      } else {
        count2=count2+1
      }
    } else {
      count1=count1+1
    }
  }
  x1 = s$level[isGeneSegment & !isAntiSenseGeneSegment]
  x2 = s$level[isGeneSegment &  isAntiSenseGeneSegment]
  cat(rt, "\n",
      "Median level of antisense genes: ", signif(median(x2), 3), "\n",
      "Median level of other genes: ", signif(median(x1), 3), "\n",
      "p=", format.pval(wilcox.test(x1, x2)$p.value), "\n\n")
}

asCand = unique(unlist(strsplit(asCand, split=", ")))

cat("The following ", length(asCand), " genes were found opposite a
segment of category ", paste("'", catgSel, "'", sep="", collapse=", "), ":\n", sep="")
print(asCand)
cat("\n\n")

cat("Analysis of 3'/5' preference of potential antisense transcripts:\n")
cat("================================================================\n")
cat("Number of cases where more than 1 feature in oppositeFeature:", count1, "\n")
cat("Number of cases where oppositeFeature not a gene:", count2, "\n")
cat("Number of cases where no unique segment for oppositeFeature:", count3, "\n")

cat("hit:\n")
print(hit)

cat("\n\nGO category analysis:\n")
cat("================================\n\n")

GOHyperG(asCand)

if(!interact) {
  dev.copy(pdf, paste(outfile, "pdf", sep="."), width=12, height=6.3)
  dev.off()
  sink()
}
