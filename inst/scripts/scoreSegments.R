##------------------------------------------------------------
## Copyright Wolfgang Huber 17 March 2005
##
## This algorithm scores segments.
## It goes sequentially through the result directory "segmentDir" in 
## ascending order, and in parallel through the GFF table.
## For each segment, we calculate and/or record its:
## * chromosome, strand, start position, end position (in bases)
## * level (mean-background)
## * t-score   (level / std error of mean)
## If it overlaps with an annotated feature
##   * the name of that feature
##   * the signed distance between its start and that of the feature
##   * the signed distance between its end and that of the feature
## If not:
##   * the distance between its start and the end of the next feature to the left
##   * the distance between its end and the start of the next feature to the left
## 
## and the same information again for the opposite strand
##------------------------------------------------------------

nrBasePerSeg = 1500 
probeLength  = 25

## For the definition of pseudogenes at SGD, see Docs/PseudogenesAtSGD.pdf
knownFeatures = c("CDS", "gene", "ncRNA", "nc_primary_transcript",
                  "rRNA", "repeat_region", "snRNA", "snoRNA", 
                  "tRNA", "transposable_element", "transposable_element_gene")

options(error=recover, warn=2)
library("tilingArray")

indir = "segmentation-050209v4"
chrs = 1:17

if(!exists("s")) {
  s  = new.env()
  totcp = 0
  cat("Loading ")
  for(chr in chrs) {
    for(strand in c("+", "-")) {
      fn = file.path(indir, paste(chr, strand, "rda", sep="."))
      cat(chr, ".", strand, " ", sep="")
      load(fn)
      assign(paste(chr, strand, "seg", sep="."), seg, envir=s)
      assign(paste(chr, strand, "dat", sep="."), dat, envir=s)
      totcp = totcp + round(max(dat$x)/nrBasePerSeg)
    }
  } ## for chr
  cat("\n")
} ## if

if(!exists("gff"))
  load("gff.rda")

pos.min = function(x) { x=x[x>=0]; if(length(x)>0) {min(x)} else {as.numeric(NA)} }

segScore = data.frame(
  chr                 = integer(totcp),
  strand              = I(character(totcp)),
  start               = integer(totcp),
  end                 = integer(totcp),
  level               = numeric(totcp),
  pt                  = numeric(totcp),
  frac.dup            = numeric(totcp),
  same.feature          = I(character(totcp)),
  same.overlap          = numeric(totcp),
  same.dist.start2feat  = integer(totcp),
  same.dist.end2feat    = integer(totcp),
  oppo.feature          = I(character(totcp)),
  oppo.overlap          = numeric(totcp),
  oppo.dist.start2feat  = integer(totcp),
  oppo.dist.end2feat    = integer(totcp))

joff = 0
for(chr in chrs) {
  for(strand in c("+", "-")) {
    dat  = get(paste(chr, strand, "dat", sep="."), s)
    seg  = get(paste(chr, strand, "seg", sep="."), s)
    cp   = round(max(dat$x)/nrBasePerSeg)
    cat(chr, ".", strand, ": ", paste(range(dat$x), collapse="..."), ", ",
        cp, " segments. ", sep="")

    ## th[i] is 1 + (end point of segment i) which is the same as
    ## the start point of segment (i-1).
    th =  c(1, seg$th[cp, 1:cp])
    i1 = th[-length(th)]         ## start points
    i2 = th[-1] - 1              ## end points

    ## idx: indices of the segments in the result table
    idx = joff + (1:cp)  
    segScore$chr[idx]      = chr
    segScore$strand[idx]   = strand
    segScore$start[idx]    = dat$x[i1]
    segScore$end[idx]      = dat$x[i2]
    segScore$frac.dup[idx] = mapply(function(h1, h2) {
      z = dat$xunique[h1:h2]
      1-sum(z)/length(z)
    }, i1, i2)

    for(wgff in c("same", "oppo")) {
      gffstrand = switch(wgff,
        same = c("+"="+", "-"="-")[strand],
        oppo = c("-"="+", "+"="-")[strand],
        stop("Sapperlot"))
      sgff = gff[ gff$seqname == chrSeqname[chr] &
                  gff$strand  == gffstrand &
                  gff$feature %in% knownFeatures, ]
      sgff$Name = getAttributeField(sgff$attributes, "Name")
      cat(wgff, ": ", nrow(sgff), "  ", sep="")

      ## matchProbes2Feats is a matrix of probes (dat$x) times features
      ## and a matrix element [p,j] is TRUE iff probe p is part of
      ## feature. A probe is considered part of a feature if its 13th
      ## nucleotide falls into it. That's an arbitrary rule, which is
      ## intended to err on the side of assigning probes to features.
      matchProbes2Feats = isWithinInterval(dat$x+(probeLength-1)/2, sgff$start, sgff$end)

      p = function(x) paste(wgff, x, sep=".")
      for(j in 1:cp) {
        js = j+joff
        ## this matrix has as many rows are there are probes in the segment,
        ## and as many columns as there as features
        matchSeg2Feats = matchProbes2Feats[i1[j]:i2[j], ]
        overlap      = colSums(matchSeg2Feats)/nrow(matchSeg2Feats)
        if(any(overlap>0)) {
          whf = which.max(overlap)
          segScore[js, p("feature")] = sgff$Name[whf]
          segScore[js, p("overlap")] = overlap[whf]
          segScore[js, p("dist.start2feat")] = segScore$start[js] - sgff$start[whf]
          segScore[js, p("dist.end2feat")]   = segScore$end[js]   - sgff$end[whf]
        } else {  
          segScore[js, p("feature")] = as.character(NA)
          segScore[js, p("overlap")] = 0
          segScore[js, p("dist.start2feat")] = pos.min(segScore$start[js] - sgff$end)
          segScore[js, p("dist.end2feat")]   = pos.min(sgff$start - segScore$end[js])
        }
      }
    }
    theCut = cut(seq(along=dat$y), th-1, labels=paste("<=", th[-1]-1, sep=""))
    
    ## use only the non-duplicated probes!
    ys     = split(dat$y[dat$xunique], theCut[dat$xunique])
    lls    = listLen(ys)
    means  = sapply(ys, mean)
    sds    = sapply(ys, sd)
    p      = pt(means/sds*sqrt(lls), df=lls-1, lower.tail=FALSE)

    ## there might not be data for all cut levels, since some have been dropped for being
    ## duplicated
    mt = match(names(means), levels(theCut))
    stopifnot(identical(names(means), names(sds)),
              identical(names(means), names(p)), !any(is.na(mt)))
    segScore$level[idx[mt]]  = means
    segScore$pt[idx[mt]]     = p

    ## just check
    ri = sample(seq(along=p), size=1)
    stopifnot(abs(p[ri] - t.test(ys[[ri]], alternative="greater")$p.value) < 1e-10)
    
    joff = joff + cp
    cat("\n")    
  }
}

save(segScore, file=file.path(indir, "segScore.rda"))
