##
## run this script after "tableSegments"
##
library("GOstats")
source("scripts/GOHyperG.R")
interact=!TRUE
if(!interact)
  sink("antisense-GO.txt")

asCand = NULL
for(rt in rnaTypes) {
  s = cs[[rt]]
  asCand = c(asCand,
    s[ s[, "category"] == "novel antisense - filtered", "oppositeFeature"])
}
asCand = unique(unlist(strsplit(asCand, split=", ")))

cat("The following", length(asCand), "genes were found opposite a
'novel antisense - filtered' segment:\n")
print(asCand)
cat("\n\n")

GOHyperG(asCand)

if(!interact) {
  dev.copy(pdf, "antisense-GO.pdf", width=12, height=6.3)
  dev.off()
  sink()
}
