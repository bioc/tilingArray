##--------------------------------------------------------------------------------
## 24.2.2005
## Use the output from the BLAST of all probes on chip vs genome to construct
## a datastructure for along-chromosome analyses
## - coord: coordinate (integer) 
##   indPM: index (1...6553600) of the PM on the chip (integer) 
##   indMM: index (1...6553600) of the MM on the chip (integer) 
##   unique: whether the probe hits multiple places (logical)
##--------------------------------------------------------------------------------
options(error=recover)

probeLength = 25
arraySize   = 2560
nrProbes    = arraySize*arraySize

if(!exists("blastRes")) {
  cat("Loading blastRes.rda.\n")
  load("blastRes.rda")
}
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

## Mapping from column 2 to chromosome ID:
## 01 = ref|NC_001133|, 02 = ref|NC_001134|, ..., 16=ref|NC_001148|
## mt = ref|NC_001224|
if(!exists("chrNo")) {
  chromoAccno = paste("ref|NC_001", c(133:148, 224), "|", sep="")
  chrNo = match(blastRes[[2]], chromoAccno)
  stopifnot(!any(is.na(chrNo)))
}

## Watson strand = antisense = "+" = typically plotted above
## Crick strand = sense = "-" = typically plotted below
if(!exists("strand")) {
  cat("Calculating cstart and cend.\n")
  cstart = blastRes[[9]]
  cend   = blastRes[[10]]
  is.antisense = (cend < cstart)    
  strand = factor(c("+", "-")[ 1 + is.antisense])
  ## swap
  tmp = cstart[is.antisense]
  cstart[is.antisense] =  cend[is.antisense] 
  cend[is.antisense] = tmp
}


if(!exists("whPM")) {
  cat("Calculating whPM: ")
  is.25mer = ((blastRes[[4]] == probeLength) &    ## alignment length
              (blastRes[[6]] == 0) &              ## number of gaps
              (blastRes[[7]] == 1) &              ## start of query
              (blastRes[[8]] == probeLength))     ## end of query

  whPM = which(is.25mer & (blastRes[[3]]==100) & (blastRes[[5]] == 0))
  whMM = which(is.25mer & (blastRes[[3]]==096) & (blastRes[[5]] == 1))
    
  ## double-check
  stopifnot(all(cend[c(whPM, whMM)]-cstart[c(whPM, whMM)]==(probeLength-1)))
  
  ## MM are always below PM!
  MMIndexFromPM = function(i) i + arraySize
  stopifnot(all(MMIndexFromPM(blastRes[[1]][whPM] %in% blastRes[[1]][whMM])))

  cat("found BLAST hits for", length(whPM), "PM probes\n")
}

## See which ones have multiple hits
PMindex = blastRes[[1]][whPM]
dupProbes = PMindex[duplicated(PMindex)]

## group PMs by chromosome and strand
sp = split(whPM, list(chr=chrNo[whPM], strand=strand[whPM]))

##--------------------------------------------------
## The along-chromosome data in 'probeAnno'
##--------------------------------------------------
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


##--------------------------------------------------
## gff
##--------------------------------------------------
if(!exists("gff")) {
  gffFile="SGD/saccharomyces_cerevisiae.gff"
  if(FALSE) {
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
    save(gff, chrSeqname, file="gff.rda", compress=TRUE)
  } else {
    load("gff.rda")
  }
}


##--------------------------------------------------
## The per-probe data in 'probeAnno'
##--------------------------------------------------
featNames = c("no_feature", "CDS", "ncRNA", "nc_primary_transcript",
              "rRNA", "snRNA", "snoRNA", "tRNA")
selGff  = which(gff$feature %in% featNames)
sgff    = gff[selGff, ]

stopifnot(nrow(sgff)>1000)
stopifnot(all(sgff$strand %in% c("+", "-")))

probe = lapply(featNames, function(i) character(nrProbes))
names(probe) = featNames


mseq = match(sgff$seqname, chrSeqname)
stopifnot(!any(is.na(mseq)), !any(mseq<1), !any(mseq>17))

for(chr in 1:length(chrSeqname)) {
  for(s in c("+", "-")) {
    ifeat  = which(sgff$strand == s & mseq == chr)
    cat("chr=", chr, "\ts=", s,
        "\tno.features=", length(ifeat), "\n", sep="")
  
    k1 = sgff$start[ifeat]
    k2 = sgff$end[ifeat]
    stopifnot(all(k1<=k2))

    p  = get(paste(chr, s, "start", sep="."), envir=probeAnno)
    k1 = matrix(k1, nrow=length(p), ncol=length(ifeat), byrow=TRUE)
    k2 = matrix(k2, nrow=length(p), ncol=length(ifeat), byrow=TRUE)

    ## sel is a logical matrix with as many rows as probes in p,
    ##   and as many columns as features in ifeat
    sel = (p>=k1 & p+(probeLength-1)<=k2)
    wh = which(colSums(sel)==0)
    if(length(wh)>0) {
      cat("\nNo probes for i=", ifeat[wh], "\n")
      print(sgff[ifeat[wh], 1:8])
    }
  
    print(table(rowSums(sel)))
 
    ind = get(paste(chr, s, "index", sep="."), envir=probeAnno)
    
    for(i in seq(along=ifeat)) {
      fi = ifeat[i]
      attrs = strsplit(sgff[fi, "attributes"], split=";")[[1]]
      g = grep("Name=", attrs)
      stopifnot(length(g)==1)
      probe[[as.character(sgff$feature[fi])]] [ind[sel[, i]]]  = sub("Name=", "", attrs[g])
    } ## for i
  } ## for s
} ## for chr


## now set 'no_feature' to all PMs that are not something else:
anyFeature = which(0 < rowSums(sapply(probe, function(s) s!="")))
allProbes = unlist(mget(paste(rep(seq(along=chrSeqname), each=2),
                              rep(c("+", "-"), length(chrSeqname)), "index", sep="."),
                   envir=probeAnno))
probe[["no_feature"]][setdiff(allProbes, anyFeature)] = "no" 
cat("length(allProbes)=", length(allProbes),
    "length(anyFeature)=", length(anyFeature),
    "no feature=", length(which(probe[["no_feature"]]=="no")), "\n")


for(i in 1:length(probe))
  probe[[i]] = factor(probe[[i]])

assign("probe", probe, envir=probeAnno)

cat("\n")



save(probeAnno, MMIndexFromPM, file="probeAnno.rda", compress=TRUE)
cat("Finished.\n")
