options(error=recover)
library("tilingArray")

indir  = "segmentation-050209v4"

if(!exists("blastres"))
  blastres = read.table(file.path(indir, "fasta", "segments.out"),
    sep="\t", as.is=TRUE, header=FALSE)

if(!exists("segScore"))
  load(file.path(indir, "segScore.rda"))
##  1  Identity of query sequence
##  2  Identity of subject sequence (matching sequence in database)
##  3  Percent identity
##  4  Alignment length
##  5  Number of mismatches
##  6  Number of gaps
##  7  Start of query sequence
##  8  End of query sequence
##  9  Start of subject sequence
## 10  End of subject sequence
## 11  E-value
## 12  Bit-score

## if(!file.exists(outdir) || !file.info(outdir)$isdir)
##   stop(paste("Output directory", outdir, "does not exist."))


## We want to distnguish three groups: annotated transcripts,
## not annotated transcripts, and not annotated not transcribed
## sequences.


groups = list(at = which(!is.na(segScore$same.feature)),
  ut = which(is.na(segScore$same.feature) & segScore$pt<1e-10),
  uu = which(is.na(segScore$same.feature) & segScore$pt>1e-2))
titgr = c("Annotated features", "Unannotated but transcribed",
          "Unannotated and not transcribed")

lE = -log(blastres[[11]], 10)
lE[lE>60] = 60
bs = blastres[[12]]
bs[bs>400] = 400

brE = seq(min(lE), max(lE), length=50)
brS = seq(min(bs), max(bs), length=50)

## A 'match' is a S. cerevisiae segment that has a BLAST hit in S. pombe.
## A 'hit' is a BLAST hit in S. pombe, corresponding to a S. cerevisiae segment
hits = lapply(groups, function(g) {
  h=which(blastres[[1]] %in% g)
  stopifnot(!any(duplicated(h)))
  return(h)
})
matches = lapply(groups, function(g) {
  which(g %in% blastres[[1]])
})

colors = brewer.pal(4, "Paired")[c(4,2,1)]
par(mfrow=c(3,2))

for(i in seq(along=groups)) {
  ntot = length(groups[[i]])
  txt = paste(length(matches), " / ", ntot, " (", signif(length(matches)/ntot, 3)*100, "%)")
  hist(lE[hits[[i]]], col=colors[i], main=titgr[i], xlab=expression(-log[10](E-value)), breaks=brE)
  hist(bs[hits[[i]]], col=colors[i], main= txt, xlab="Bit-score", breaks=brS)
}

dev.copy(pdf, file="segConservation.pdf", width=7, height=10); dev.off()

print(wilcox.test(bs[hits[[2]]], bs[hits[[1]]]))
print(wilcox.test(bs[hits[[2]]], bs[hits[[3]]]))

