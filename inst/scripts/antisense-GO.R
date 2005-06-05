

##
## run this script after "tableSegments"
##

library("GO")
library("annotate")

interact=!TRUE
if(!interact)
  sink("antisense-GO.txt")

if(!("Ontology_term" %in% names(gff)))
  gff$Ontology_term=getAttributeField(gff[, "attributes"], "Ontology_term")

asCand = NULL
for(rt in rnaTypes) {
  s = cs[[rt]]
  asCand = c(asCand,
    s[ s[, "category"] == "novel antisense - filtered", "oppositeFeature"])
}
asCand = unique(unlist(strsplit(asCand, split=", ")))

ctrls  = gff[ gff[, "feature"]=="gene", "Name"]
ctrls  = unlist(strsplit(ctrls, split=", "))
ctrls  = setdiff(ctrls, asCand)

cat("The following", length(asCand), "genes were found opposite a
'novel antisense - filtered' segment:\n")
print(asCand)
cat("\n\n")

## get the GO classes for a list of genes:
whg = which(gff[, "feature"]=="gene")
getGO = function(x) {
  mt  = match(x, gff[whg, "Name"])
  #if(any(is.na(mt)))
  #  cat("The following have no match:", x[is.na(mt)], "\n")
  rv = strsplit(gff[whg[mt], "Ontology_term"], split=",")
  stopifnot(!any(sapply(rv, function(x) any(duplicated(x)))))
  rv
}

asGO = getGO(asCand)
ctGO = getGO(ctrls)

asTab = table(unlist(asGO))
ctTab = table(unlist(ctGO))

allGO = sort(unique(c(unlist(asGO), unlist(ctGO))))
nr = matrix(NA, nrow=length(allGO), ncol=2)
rownames(nr) = allGO

nr[names(asTab), 1] = asTab
nr[names(ctTab), 2] = ctTab

plot(nr, pch=16, xlab="antisense genes", ylab="control: all genes",
     main="Frequency of occurence of GO classes")
thSlop = length(ctrls) / length(asCand)
abline(a=0, b=thSlop, col="blue")
dev.copy(pdf, "antisense-GO.pdf", width=7, height=7); dev.off()


for(k in c(48, 9, 7:5)) {
  goL = rownames(nr)[ which(nr[,1]==k & nr[,2]/nr[,1]<thSlop) ]

  cat("GO classes which occured ", k, " times (and below the blue line).\n",
      "============================================================\n", sep="")
  for(L in goL) {
    print(sort(replaceSystematicByCommonName(asCand[sapply(asGO, function(x) L %in% x)])))
    print(get(L, GOTERM))
    cat("\n")
  }
}

if(!interact)
  sink()
