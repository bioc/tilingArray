##
## How many known features do we find exprssed?
##
library("tilingArray")
library("multtest")

rnaTypes  = c("seg-polyA-050525", "seg-tot-050525")
source(scriptsDir("readSegments.R"))

interact=!TRUE
if(!interact) {
  sink("expressedFeatures.txt")
  cat("Made on", date(), "\n\n")
}

source(scriptsDir("calcThreshold.R"))

if(!exists("intergenic")){
  probe = probeAnno$probeReverse
  intergenic = which(probeAnno$probeReverse$no_feature=="no" & probeAnno$probeDirect$no_feature=="no")
}

for(rt in rnaTypes) {
  e = get(rt)
  if(!("xn" %in% ls(e))) {
    fn = file.path(rt, "xn.rda")
    cat("Loading", fn, "\n")
    load(fn, envir=e)
    x = exprs(get("xn", e))
    for(j in 1:ncol(x))
      x[,j] = x[, j] - median(x[intergenic, j])
    assign("x", x, envir=e)
    rm(x)
  }
}

sel = gff[, "feature"]=="gene"
allncRNA = c("ncRNA","snoRNA","snRNA","tRNA","rRNA")

if(!"ID"%in%colnames(gff))
  gff$ID=getAttributeField(gff[, "attributes"], "ID")

sgff = rbind(
  cbind(category="verified gene",        gff[ sel & gff[, "orf_classification"]=="Verified", ]),
  cbind(category="uncharacterized gene", gff[ sel & gff[, "orf_classification"]=="Uncharacterized", ]),
  cbind(category="dubious gene",         gff[ sel & gff[, "orf_classification"]=="Dubious", ]),
  cbind(category="ncRNA(all)",           gff[ gff[, "feature"] %in% allncRNA & !is.na(gff[, "ID"]), ]))
colnames(sgff)[1]="category"

cat("Features in GFF table:\n")
print(table(sgff$category))

## duplicated?
multip = split(1:nrow(sgff), sgff[, "Name"])
for(w in which(listLen(multip)>1))
  for(cn in c("chr", "start", "end", "strand"))
    stopifnot(length(unique(sgff[multip[[w]], cn]))==1)

stopifnot(!any(is.na(sgff[, "Name"])), !any(duplicated(sgff[, "Name"])))

## for each interesting feature (gene, RNA), build a list of probes:
if(!exists("index")) {
  index = vector(mode="list", length=nrow(sgff))
  names(index) = sgff[, "Name"]
  
  ind = paste(gff$chr, gff$strand, "index", sep=".")
  sta = paste(gff$chr, gff$strand, "start", sep=".")
  end = paste(gff$chr, gff$strand, "end",   sep=".")
  uni = paste(gff$chr, gff$strand, "unique",   sep=".")
  stopifnot(is.numeric(uni))
  
  sp = split(1:nrow(gff), gff[, "Name"])
  for(i in seq(along=index)) {
    ## this combination of for/if is really slow
    whf = sp[[names(index)[i]]]
    fn  = gff[whf, "feature"]
    ## genes (may contain introns)
    if(sum(fn=="gene")==1 && all(fn%in%c("gene", "CDS", "intron", "region"))) {
      fset = which(fn=="CDS")
      ## tRNAs (may contain introns)  
    } else if (sum(fn=="tRNA")==1 && all(fn%in%c("tRNA","ncRNA","intron"))) {
      fset = which(fn=="ncRNA")
      ## rRNAs (may contain introns)  
    } else if (sum(fn=="rRNA")==1 && all(fn%in%c("rRNA","ncRNA","intron", "nc_primary_transcript"))) {
      fset = which(fn%in%c("ncRNA","nc_primary_transcript"))
      ## snoRNAs (may contain introns)  
    } else if (sum(fn=="snoRNA")==1 && all(fn%in%c("snoRNA","ncRNA", "intron"))) {
      fset = which(fn=="ncRNA")
      ## snRNAs   
    } else if (sum(fn=="snRNA")==1 && all(fn%in%c("snRNA","ncRNA"))) {
      fset = which(fn=="ncRNA")
    } else if (all(fn=="ncRNA")) {
      fset = which(fn=="ncRNA")
    } else {
      cat("Dropping:\n")
      print(gff[whf, c(1, 3:5, 6)])
      cat(gsub("%20", " ", unique(getAttributeField(gff[whf, "attributes"], "Note"))), "\n\n")
      fset = integer(0)
    }
    ## gene with one or more CDSs
    res = integer(0)
    for(w in whf[fset])
      res = c(res, get(ind[w], probeAnno)[
        (get(uni[w], probeAnno)==0) &
        (get(sta[w], probeAnno) >= gff$start[w])&
        (get(end[w], probeAnno) <= gff$end[w]) ])
    index[[i]] = sort(unique(res))
  }
}

hasEnoughProbes = (listLen(index)>=7)
cat(sum(hasEnoughProbes),"of",length(index),"potential transcripts are matched by >=7 unique probes.\n\n")

res = matrix(NA, nrow=length(levels(sgff[,"category"])), ncol=8)
stopifnot(all(rnaTypes==c("seg-polyA-050525", "seg-tot-050525")))
rownames(res) = levels(sgff[,"category"])
colnames(res) = c("n1: in genome", "n2: with probes", 
    "det. poly-A", "poly-A: % of n1", "poly-A: % of n2",
    "det. total",  "total: % of n1",  "total: % of n2")

res[names(n1), 1] = n1 = table(sgff[, "category"])
res[names(n2), 2] = n2 = table(sgff[hasEnoughProbes, "category"])


for(irt in seq(along=rnaTypes)) {
  rt = rnaTypes[irt]
  ## cat("\n=====", rt, "=====\n")
  x = get(rt)$"x"

  doTest = function(ind) {
    pv = sapply(ind, function(jj)
           binom.test(sum(x[jj, ]>0), ncol(x)*length(jj), alternative="greater")$p.value)
    
    bh = mt.rawp2adjp(pv, proc="BY")
    stopifnot(all(bh$adjp[, 1] == pv[bh$index], na.rm=TRUE))
    adjp = numeric(nrow(bh$adjp))
    adjp[bh$index] = bh$adjp[,2]
    (adjp < FDRthresh)
  }

  isExp = doTest(index[hasEnoughProbes])
      
  nrDet = table(sgff[hasEnoughProbes, "category"], isExp)[, "TRUE"]

  stopifnot(identical(names(nrDet), names(n1)), identical(names(nrDet), names(n2)))
  stopifnot(setequal(names(nrDet), rownames(res)))
  
  res[names(nrDet), irt*3]   = nrDet
  res[names(nrDet), irt*3+1] = round(100*nrDet/n1, 1)
  res[names(nrDet), irt*3+2] = round(100*nrDet/n2, 1)    
}
  
print(res)

##
##  GO analysis of unexpressed genes
##
cat("\n\nGO-Analysis of the untranscribed verified genes:\n\n")

source(scriptsDir("GOHyperG.R"))
source(scriptsDir("writeSegmentTable.R"))

unexpressedGenes = names(index)[hasEnoughProbes][!isExp & sgff[hasEnoughProbes, "category"]=="verified gene"]


GOHyperG(unexpressedGenes)

if(!interact)
  sink()
