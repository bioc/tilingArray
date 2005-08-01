##
## run this script after "tableSegments"
##
library("GOstats")
source("scripts/GOHyperG.R")

what = c("filtered", "all")[2]
switch(what,
       filtered = {
         outfile = "antisense-GO"
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
for(rt in rnaTypes) {
  s = cs[[rt]]
  asCand = c(asCand,
    s[ s[, "category"] %in% catgSel, "oppositeFeature"])
}
asCand = unique(unlist(strsplit(asCand, split=", ")))

cat("The following ", length(asCand), " genes were found opposite a
segment of category ", paste("'", catgSel, "'", sep="", collapse=", "), ":\n", sep="")
print(asCand)
cat("\n\n")

GOHyperG(asCand)

if(!interact) {
  dev.copy(pdf, paste(outfile, "pdf", sep="."), width=12, height=6.3)
  dev.off()
  sink()
}
