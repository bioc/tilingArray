##--------------------------------------------------------------------------------
## Use the output from the BLAST of all probes on the chip vs genome to construct
## probe annotaton datastructures
##
## 26.3.2005:
## Note: this script is such a memory hog that it currently only works on
## dijkstra (27 GB). On pavlov (16 GB) it will fail.
##
## 13.5.2005:
## Extended to also write arrayDesign.txt file for ArrayExpress
##
## The script has the following parts
## 1. Read and process blastRes
## 2. Write arrayDesign.txt
## 3. Construct along-chromosome vectors for probeAnno environment:
##    coord: coordinate (integer) 
##    index: index (1...6553600) of the PM on the chip (integer) 
##    unique: whether the probe hits multiple places (logical)
## 4. Construct GFF dataframe
## 5. Construct probeAnnoReverse and probeAnnoDirect lists
## 6. Write 3-5 into probeAnno.rda
##--------------------------------------------------------------------------------
library("tilingArray")
library("Scerevisiaetilingprobe")

options(error=recover)

probeLength = 25
arraySize   = 2560
nrProbes    = arraySize*arraySize

##
## PART 1
##
## blastRes is a dataframe with  13,550,085 rows and 12 columns
##  1  Identity of query sequence  (integer)
##  2  Identity of subject sequence = matching sequence in database (character)
##  3  Percent identity (numeric)
##  4  Alignment length (integer)
##  5  Number of mismatches (integer)
##  6  Number of gaps (integer)
##  7  Start of query sequence (integer)
##  8  End of query sequence (integer)
##  9  Start of subject sequence (integer)
## 10  End of subject sequence (integer)
## 11  E-value (numeric)
## 12  Bit-score (numeric)
if(!exists("blastRes")) {
  cat("Loading blastRes.rda.\n")
  load("blastRes.rda")
}

## Mapping from column 2 to chromosome ID:
## 01 = ref|NC_001133|, 02 = ref|NC_001134|, ..., 16=ref|NC_001148|
## mt = ref|NC_001224|
if(!exists("chrNo")) {
  chromoAccno = paste("ref|NC_001", c(133:148, 224), "|", sep="")
  chrNo = match(blastRes[[2]], chromoAccno)
  stopifnot(!any(is.na(chrNo)))
}

## Watson strand = "+" = typically plotted above
## Crick strand  = "-" = typically plotted below
if(!exists("strand")) {
  cat("Calculating cstart and cend.\n")
  cstart = blastRes[[9]]
  cend   = blastRes[[10]]
  isMinusStrand = (cend < cstart)    
  strand = factor(c("+", "-")[ 1 + isMinusStrand])
  ## swap
  tmp = cstart[isMinusStrand]
  cstart[isMinusStrand] =  cend[isMinusStrand] 
  cend[isMinusStrand] = tmp
}


if(!exists("whPM")) {
  cat("Calculating whPM: ")
  is.25mer = ((blastRes[[4]] == probeLength) &    ## alignment length
              (blastRes[[6]] == 0) &              ## number of gaps
              (blastRes[[7]] == 1) &              ## start of query
              (blastRes[[8]] == probeLength))     ## end of query

  whPM = which(is.25mer & (blastRes[[3]]==100) & (blastRes[[5]] == 0))
  whMM = which(is.25mer & (blastRes[[3]]==096) & (blastRes[[5]] == 1))


  stop("INTERRUPT")
  
  ## double-check
  stopifnot(all(cend[c(whPM, whMM)]-cstart[c(whPM, whMM)]==(probeLength-1)))
  
  ## MM are always below PM!
  MMIndexFromPM = function(i) i + arraySize
  stopifnot(all(MMIndexFromPM(blastRes[[1]][whPM] %in% blastRes[[1]][whMM])))

  ## See which ones have multiple hits
  PMindex = blastRes[[1]][whPM]
  dupProbes = PMindex[duplicated(PMindex)]
  
  ## group PMs by chromosome and strand
  sp = split(whPM, list(chr=chrNo[whPM], strand=strand[whPM]))
  
  cat("found", length(whPM), "perfect match and", length(whMM),
      "mismatch BLAST hits.\n")
}

##
## Part 2
## 
  
idx = as.integer(Scerevisiaetilingprobe$y*2560 + Scerevisiaetilingprobe$x + 1)
stopifnot(identical(idx, 1:length(Scerevisiaetilingprobe$x)))
            
mt = whPM[match(idx, blastRes[[1]][whPM])]

out = c(paste("Row","Column","Reporter.name",
  "Reporter.Sequence",
  "Target.Name",
  "Target.Start",
  "Target.End",
       sep="\t"), 
  paste(Scerevisiaetilingprobe$y,
        Scerevisiaetilingprobe$x,
        idx,
        reverseSeq(Scerevisiaetilingprobe$sequence),
        blastRes[[2]][mt],
        blastRes[[9]][mt],
        blastRes[[10]][mt],
  sep="\t"))

con = gzfile("ForArrayExpress/arrayDesign.txt.gz", open="wt")
writeLines(out, con=con)
close(con)

stop("Stopped after Part 2.")

##--------------------------------------------------
##
## Part 3
##
## The along-chromosome data in 'probeAnno'
##--------------------------------------------------
if(!exists("probeAnno")) {
  cat("Making probeanno: ")
  probeAnno = new.env()
  for(i in seq(along=sp)) {
    whc = sp[[i]]
    nm  = names(sp)[i]
    cat(nm, "")
    PMindex = blastRes[[1]][whc]
    assign(paste(nm, "start", sep="."),  cstart[whc], envir=probeAnno)
    assign(paste(nm, "index", sep="."),  PMindex, envir=probeAnno)
    assign(paste(nm, "unique", sep="."), !(PMindex %in% dupProbes), envir=probeAnno)
  }
  cat("\n")
}

##--------------------------------------------------
## Part 4: gff
##--------------------------------------------------
if(!exists("gff")) {
  gffFile="SGD/saccharomyces_cerevisiae.gff"
  cat("Reading gff: ")
  ## GFF Files
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="", header=FALSE, comment.char="#",
    colClasses=c("factor", "factor", "factor", "integer", "integer",
      "factor", "factor", "factor", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
            "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes", paste(sapply(gff, class), collapse=" "), ".\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  
  chrSeqname = 
    c("chrI", "chrII", "chrIII", "chrIV",
      "chrV", "chrVI", "chrVII", "chrVIII", "chrIX",
      "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV",
      "chrXV","chrXVI", "chrMito")
}

##--------------------------------------------------
## Part 5: the per-probe data in 'probeAnno'
##--------------------------------------------------
featNames = c("no_feature", "CDS", "ncRNA", "nc_primary_transcript",
              "rRNA", "snRNA", "snoRNA", "tRNA")
selGff  = which(gff$feature %in% featNames)  
sgff    = gff[selGff, ]

stopifnot(nrow(sgff)>1000)
stopifnot(all(sgff$strand %in% c("+", "-")))  

mseq = match(sgff$seqname, chrSeqname)
stopifnot(!any(is.na(mseq)), !any(mseq<1), !any(mseq>17))

for(hybeType in c("Reverse", "Direct")) {
  cat("-----", hybeType, "-----\n")
  p = lapply(featNames, function(i) character(nrProbes))
  names(p) = featNames
  ## loop over chromosomes and strands
  for(chr in 1:length(chrSeqname)) {
    for(s in c("+", "-")) {
      ifeat  = which(sgff$strand == s & mseq == chr)
      cat(chr, s, ": ", length(ifeat), " features\n", sep="")

      k1 = sgff$start[ifeat]
      k2 = sgff$end[ifeat]
      stopifnot(all(k1<=k2))
   
      Name = getAttributeField(sgff$attributes[ifeat], "Name")

      pas = switch(hybeType,
        Reverse = s,
        Direct  = otherStrand(s)
      )
      sta = get(paste(chr, pas, "start", sep="."), envir=probeAnno)
      ind = get(paste(chr, pas, "index", sep="."), envir=probeAnno)
      sel = isWithinInterval(sta+(probeLength-1)/2, k1, k2)
      wh  = which(colSums(sel)==0)
      ## if(length(wh)>0) {
      ##	cat("\nNo probes for i=", ifeat[wh], "\n")
      ##	print(sgff[ifeat[wh], 1:8])
      ## }
      ## print(table(rowSums(sel)))
       
      for(i in seq(along=ifeat)) {
	fi = ifeat[i]
	p[[as.character(sgff$feature[fi])]] [ind[sel[, i]]]  = Name[i]
      } ## for i
    } ## for s
  } ## for chr

  ## now set 'no_feature' to all PMs that are not something else:
  anyFeature = which(0 < rowSums(sapply(p, function(s) s!="")))
  allProbes = unlist(mget(paste(rep(seq(along=chrSeqname), each=2),
                              rep(c("+", "-"), length(chrSeqname)), "index", sep="."),
    envir=probeAnno))
  p[["no_feature"]][setdiff(allProbes, anyFeature)] = "no" 
  cat("length(allProbes)=", length(allProbes),
      "length(anyFeature)=", length(anyFeature),
      "no feature=", length(which(p[["no_feature"]]=="no")), "\n")

  for(i in 1:length(p))
    p[[i]] = factor(p[[i]])

  assign(paste("probe", hybeType, sep=""), p, envir=probeAnno)
} ## for hybetype

cat("\n")
##
##  Part 6
##
save(probeAnno, MMIndexFromPM, gff, chrSeqname, file="probeAnno.rda", compress=TRUE)
cat("Finished.\n")
