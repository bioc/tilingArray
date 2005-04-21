# compare signal to noise ratios before and after DNA normalization

library("tilingArray")
library("genefilter")
library("multtest")
pthresh=0.05

if(!exists("a"))load("a.rda")
if(!exists("x"))load("x.rda")
if(!exists("gff")) {
  load("probeAnno.rda")
  gff$Name = getAttributeField(gff$attributes, "Name")
  gff$orf_classification = getAttributeField(gff$attributes, "orf_classification")

  ## the names of the genes in three classes:  Dubious, Uncharacterized, Verified
  orfClasses = unique(gff$orf_classification)
  orfClasses = orfClasses[!is.na(orfClasses)]
  geneNames = lapply(orfClasses, function(cl)
    gff$Name[gff$feature=="gene" & gff$orf_classification==cl])
  names(geneNames) = orfClasses
  
  ## standard deviations for annotated CDSs
  probe=probeAnno$probeReverse
}

## x can be a matrix
calcSds = function(x) {
  stopifnot(is.matrix(x))
  for(i in 1:ncol(x)) {
    cat(i, "")
    xc = x[,i]
    ## baseline = shorth(xc[probe$no_feature=="no"])
    baseline = median(xc[probe$no_feature=="no"])
    xc = xc-baseline

    xs  = split(xc, probe$CDS)
    if(i==1) {
      ## throw out CDS with less than 7 probes
      xs  = xs[(names(xs)!="" & listLen(xs)>=7)]
      means = sds = matrix(as.numeric(NA), ncol=ncol(x), nrow=length(xs))
      rownames(means) = rownames(sds) = names(xs)
    } else {
      xs = xs[rownames(means)]
    }
    
    means[,i] = sapply(xs, mean)
    sds[,i]   = sapply(xs, sd)
  } ## for i

  nr    = ncol(sds)*listLen(xs)
  p     = pt( rowMeans(means) / sqrt(rowSums(sds*sds)) * sqrt(nr), df=nr-1, lower.tail=FALSE)
  bh    = mt.rawp2adjp(p, proc="BY")
  stopifnot(all(bh$adjp[, 1] == p[bh$index]))
  
  p[bh$index] = bh$adjp[,2]
  names(p) = rownames(means)
  
  list(mean=means, sd=sds, p=p)
}

fn = list("polyA" = 
  c("041203_S96_polyAx1_RH6.cel.gz",
    "050209_mRNAx4_30min_re-hybe_RH6.cel.gz",
    "050218_polyA-RNA_RH6_4x15min.cel.gz"),
  "tot" = 
  c("050409_totcDNA_14ug_no52.cel.gz",
    "030505_totcDNA_15ug_affy.cel.gz"))
    ##"050415_totcDNA_20ug_Affy11.cel.gz"))

if(!exists("res")) 
  res = lapply(fn, function(f) {
    before = calcSds(log(exprs(a)[, f], 2))
    after  = calcSds(exprs(x)[, f])
    list(before=before, after=after)
  })


out = stdout()
cat("Normalization:\n", file=out)

par(mfrow=c(2,1))
from=0; to=3; n=512; cols=c("black", "blue")

for(rt in seq(along=fn)) {
  nm = c("poly-A RNA", "total RNA")[rt]

  ## plot the standard deviations before and after
  mean.before  = res[[rt]]$before$mean
  sel          = (mean.before>median(mean.before))
  before    = res[[rt]]$before$sd[sel]
  after     = res[[rt]]$after$sd[sel]
  
  db = density(before, from=from, to=to, n=n)
  n  = length(db$x)
  da = density(after, from=from, to=to, n=n)
  ymax = max(db$y, da$y)
  plot(db, ylim=c(0,ymax), lwd=2, col=cols[1], main=nm, xlab="standard deviation")
  lines(da, lwd=2, col=cols[2])
  if(rt==1)
    legend(2.2, 0.8*ymax, legend=c("before", "after"), text.width=0.5, y.intersp=2, lty=1, col=cols, lwd=2)

  cat("\n\n", nm, ": median sds before=", signif(median(before), 3),
      ", after=", signif(median(after), 3),
      ", median reduction=", signif(median(before/after), 3), "\n", file=out, sep="")

  ## print the number of detected genes
  p = res[[rt]]$after$p
  ## stopifnot(all(names(p) %in% unlist(geneNames)))
  for(cl in orfClasses) {
    mt  = match(geneNames[[cl]], names(p))
    pcl = p[mt[!is.na(mt)]]
    ntot = length(mt)
    ndet = sum(pcl<=0.05)
    cat("ORF class '", cl, "' contains ", ntot, " genes of which we have p-values for ",
        sum(!is.na(mt)), " and ", ndet, " are <=", pthresh, " (", signif(ndet/ntot*100, 3), 
        "%).\n", sep="", file=out)
  }
}


dev.copy(pdf, file="fig_norm.pdf", width=6, height=7); dev.off()

