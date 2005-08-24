##--------------------------------------------------------------------------------
## Use the output from the MUMmer (querying all probes on the chip vs genome)
## to construct probe annotaton datastructures
##
## 06.8.2005: Use MUMmer instead of BLAST
## 21.5.2005: Use also probes shorter than 25 bases
## 13.5.2005: Write arrayDesign.txt file for ArrayExpress
##
## The script has the following parts
## 1. Read and process scAll.out
## 2. Write arrayDesign.txt
## 3. Construct along-chromosome vectors for probeAnno environment:
##    start, end: (integer) 
##    index: index (1...6553600) of the PM on the chip (integer) 
##    unique: whether the probe hits multiple places (logical)
## 4. Construct GFF dataframe
## 5. Construct probeAnnoReverse and probeAnnoDirect lists
## 6. Write 3-5 into probeAnno.rda
##--------------------------------------------------------------------------------
library("tilingArray")

options(error=recover)

arraySize   = 2560
nrProbes    = arraySize*arraySize

what = c("hits", "arrayexpress", "probeAnno1", "gff", "probeAnno2", "probeAnnoSave")[5:6]

if("hits" %in% what) {
  ## read the MUMmer output, parse and process
  mmRes = readLines("SGD-0508/scAll.out")

  ## ------------------------------------------------------------------------
  ## parse the sequence header lines and read them into three integer vectors
  ## seqId, seqDir, seqLen
  ## ------------------------------------------------------------------------
  seqj  = grep("^>", mmRes) 
  
  seqLines = sub("Reverse  Len =", "-1", mmRes[seqj])
  seqLines = sub(" Len =", "1", seqLines)
  seqLines = sub("> ", "", seqLines)
  sp  = strsplit(seqLines, " ")
  stopifnot(all(listLen(sp)==3))
  spi = matrix(as.integer(unlist(sp)), nrow=3)
  seqLen   = spi[3,]
  seqDir   = spi[2,]
  seqId    = spi[1,]
  
  odd  = seq(1, length(seqId), by=2)
  even = seq(2, length(seqId), by=2)
  stopifnot(all(seqDir==c(1,-1)), all(seqId[odd]==seqId[even]))
  
  rm(list=c("seqLines", "seqDir", "sp", "spi"))
  
  chromoAccno = paste("ref|NC_001", c(133:148, 224), "|", sep="")
  
  ## --------------------------------------------------------------------
  ## parse the hit lines and put results into matrix 'hit', with columns:
  ##  seqId:      1..65536000
  ##  chr:        1..17
  ##  strand:     -1, +1 (+1 for Watson/sense, -1 for Crick/antisense)
  ##  start:      1...n. "start" is always the numerically lowest
  ##              coordinate, i.e. the 5' base for Watson strand, the 3'
  ##              base for Crick strand!
  ##  length:     1...L
  ##  unique:     3: has multiple perfect matches (PMs) in the genome
  ##              2: has no PM but one or more near-matches 
  ##              1: has exactly one PM and some near-matches in the genome
  ##              0: has exactly one PM and no near-matches
  ## Here, what is a near match is determined by the parameter settings of
  ## mummer, see 'mapProbesToGenome.sh'. In particular, with he parameter
  ## "-l 23" we require at least a perfect match of 23 probes
  ## ---------------------------------------------------------------------
  
  hits = matrix(as.integer(NA), nrow=length(mmRes)-length(seqLen), ncol=6)
  colnames(hits) = c("seqId", "chr", "strand", "start", "length", "unique")
  k = 0
  jstart = seqj+1
  jend   = c(seqj[-1]-1, length(mmRes)-1)
  
  ## this for-loop will run for about 1-2h...
  for(i in odd) {
    if(i%%10000==1)
      cat(i, "")
    jf = jstart[i]  : jend[i]
    jr = jstart[i+1]: jend[i+1]
    if(jf[length(jf)]<jf[1])
      jf=integer(0)
    if(jr[length(jr)]<jr[1])
      jr=integer(0)
    j  = c(jf, jr)
    if(length(j)>0) {
      sl      = seqLen[i]
      strand  = rep(c(+1, -1), c(length(jf), length(jr)))
      sp      = strsplit(mmRes[j], " +", perl=TRUE)
      stopifnot(all(listLen(sp)==5), all(sapply(sp, "[", 1)==""))
      chr     = match(sapply(sp, "[", 2), chromoAccno)
      stopifnot(!any(is.na(chr)))
      start   = as.integer(sapply(sp, "[", 3))
      posQur  = as.integer(sapply(sp, "[", 4))
      lenMt   = as.integer(sapply(sp, "[", 5))
      isPM    = (lenMt==sl)
      stopifnot( all (!isPM | ((sl/2+0.5)-(sl/2-0.5)*strand==posQur)))
      uniq    = if(sum(isPM)==1) {
                  as.integer(!all(isPM))
                } else {
                  2 + (sum(isPM)>1)
                }
      hits[k+seq(along=j), ] = cbind(seqId[i], chr, strand, start, sl, uniq)
      k = k+length(j)
    }
  }
  stopifnot(k==nrow(hits))
  save(hits, file="makeProbeAnno-hits.rda")
}  else {
  ## read the cached file "hits"
  if(!exists("hits")) {
    cat("Loading makeProbeAnno-hits.rda\n")
    load("makeProbeAnno-hits.rda")
  }
} ## if 

if(!exists("whPM")) {
  whPM = which(hits[, "unique"]!=2)
  print(table(hits[,"unique"]))
  ##       0       1       2       3 
  ## 2834296   42118 2099819  937386 
}


##
## Part 2
##
if("arrayexpress" %in% what) {
  library("Scerevisiaetilingprobe")
  aefn = "ForArrayExpress/arrayDesign-0508.txt.gz"
  cat("Writing", aefn, "\n")
  idx = as.integer(Scerevisiaetilingprobe$y*2560 + Scerevisiaetilingprobe$x + 1)
  stopifnot(identical(idx, 1:length(Scerevisiaetilingprobe$x)))
  
  mt  = whPM[match(idx, hits[whPM, "seqId"])]
  out = c(paste("Row","Column","Reporter.name",
    "Reporter.Sequence",
    "Target.Name",
    "Target.Start",
    "Target.Strand",
    "Reporter.Status",
    sep="\t"), 
    paste(Scerevisiaetilingprobe$y,
          Scerevisiaetilingprobe$x,
          idx,
          ifelse(is.na(Scerevisiaetilingprobe$sequence), "",
                 reverseSeq(Scerevisiaetilingprobe$sequence)),
          chromoAccno[ hits[mt, "chr"] ],
          hits[mt, "start"],
          hits[mt, "length"],
          hits[mt, "unique"],
          sep="\t"))
  
  con = gzfile(aefn, open="wt")
  writeLines(out, con=con)
  close(con)

}


##--------------------------------------------------
##
## Part 3
##
## The along-chromosome data in 'probeAnno'
##--------------------------------------------------
if("probeAnno1" %in% what) { 
  ## group PMs by chromosome and strand
  cat("Making probeanno: ")
  sp = split(whPM, list(chr=hits[whPM, "chr"], strand=c("-", "+")[(3+hits[whPM, "strand"])/2]))
  probeAnno = new.env()
  for(i in seq(along=sp)) {
    whc = sp[[i]]
    nm  = names(sp)[i]
    cat(nm, "")
    assign(paste(nm, "start", sep="."),  as.integer(hits[whc, "start"]), envir=probeAnno)
    assign(paste(nm, "end", sep="."),    as.integer(hits[whc, "start"]+hits[whc, "length"]-1),
            envir=probeAnno)
    assign(paste(nm, "index", sep="."),  as.integer(hits[whc, "seqId"]),  envir=probeAnno)
    assign(paste(nm, "unique", sep="."), as.integer(hits[whc, "unique"]), envir=probeAnno)
  }
  cat("\n")
}

##--------------------------------------------------
## Part 4: gff
##--------------------------------------------------
if("gff" %in% what) { 
  nrchr = 17
  ## GFF Files
  gffRead = function(gffFile) {
    cat("Reading ", gffFile, ": ", sep="")
    gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="", header=FALSE, comment.char="#",
      colClasses=c("factor", "factor", "factor", "integer", "integer",
        "factor", "factor", "factor", "character"))
    colnames(gff) = c("seqname", "source", "feature", "start", "end",
              "score", "strand", "frame", "attributes")
    cat("found", nrow(gff), "rows with classes:", paste(sapply(gff, class), collapse=", "), "\n")
    stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
    return(gff)
  }
  
  ## also read the regulatory features:
  gff1 = gffRead("SGD-0508/saccharomyces_cerevisiae.gff")
  gff2 = gffRead("SGD-0508/scerevisiae_regulatory.gff")
  gff  = rbind(gff1, gff2)

  ## Add additional useful fields
  gff$Name = getAttributeField(gff$attributes, "Name")
  theID    = getAttributeField(gff$attributes, "ID")
  stopifnot(all(gff$Name == theID, na.rm=TRUE))
  gff$Name[is.na(gff$Name)] <- theID[is.na(gff$Name)]
  
  gff$orf_classification = getAttributeField(gff$attributes, "orf_classification")
  gff$gene               = getAttributeField(gff$attributes, "gene")

  gff$chr = match(gff$seqname,
    c("chrI", "chrII", "chrIII", "chrIV",
      "chrV", "chrVI", "chrVII", "chrVIII", "chrIX",
      "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV",
      "chrXV","chrXVI", "chrMito"))
  stopifnot(!any(is.na(gff$chr)), !any(gff$chr<1), !any(gff$chr>nrchr))

  # 2005/08/24: added "CDS_dubious" as new (more precise) feature description
  gff$feature <- as.character(gff$feature)
  gff$feature[gff$feature=="CDS" & gff[,"orf_classification"]=="Dubious"] <- "CDS_dubious"

}

##--------------------------------------------------
## Part 5: the per-probe data in 'probeAnno'
##--------------------------------------------------
if("probeAnno2" %in% what) { 
featNames = c("no_feature", "CDS", "ncRNA", "nc_primary_transcript",
              "rRNA", "snRNA", "snoRNA", "tRNA")
selGff  = which(gff$feature %in% featNames)  
sgff    = gff[selGff, ]

stopifnot(nrow(sgff)>1000)
stopifnot(all(sgff$strand %in% c("+", "-")))  

for(hybeType in c("Reverse", "Direct")) {
  cat("-----", hybeType, "-----\n")
  p = lapply(featNames, function(i) character(nrProbes))
  names(p) = featNames

  ## loop over chromosomes and strands
  for(s in c("+", "-")) {
    pas = switch(hybeType,
      Reverse = s,
      Direct  = otherStrand(s)
      )
    for(chr in 1:nrchr) {
      ifeat  = which(sgff$strand == s & sgff$chr == chr)
      cat(chr, s, ": ", length(ifeat), " features\n", sep="")

      k1   = sgff[ifeat, "start"]
      k2   = sgff[ifeat, "end"]
      Name = sgff[ifeat, "Name"]
      feat = as.character(sgff[ifeat, "feature"])
      stopifnot(all(k1<=k2))

      sta = get(paste(chr, pas, "start", sep="."), envir=probeAnno)
      end = get(paste(chr, pas, "end",   sep="."), envir=probeAnno)
      ind = get(paste(chr, pas, "index", sep="."), envir=probeAnno)

      for(i in seq(along=ifeat)) {
        sel = (sta>=k1[i] & end<=k2[i])
	p[[feat[i]]][ind[sel]] = Name[i]
      } ## for i
    } ## for chr
  } ## for s

  ## now set 'no_feature' to all PMs that are not something else:
  anyFeature = which(0 < rowSums(sapply(p, function(s) s!="")))
  allProbes = unique(unlist(mget(paste(
    rep(1:nrchr, each=2), rep(c("+", "-"), nrchr), "index", sep="."),
    envir=probeAnno)))
  
  p[["no_feature"]][setdiff(allProbes, anyFeature)] = "no" 
  cat("length(allProbes)=", length(allProbes),
      "length(anyFeature)=", length(anyFeature),
      "no feature=", length(which(p[["no_feature"]]=="no")), "\n")

  for(i in 1:length(p))
    p[[i]] = factor(p[[i]])

  assign(paste("probe", hybeType, sep=""), p, envir=probeAnno)
} ## for hybetype
cat("\n")
}

##
##  save
##
if("probeAnnoSave" %in% what) { 
  cat("Saving probeAnno.rda.\n")
  save(probeAnno, gff, file="probeAnno.rda", compress=TRUE)
}
