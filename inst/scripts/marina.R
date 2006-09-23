## A script for Marina Granovskaia to make Along-Chromosome Plots
## (C) W. Huber 2006
##
## Please specify the parameters in the following lines:

## Name of RNA and of DNA hybes.
## RNA must be a single files, for DNA multiple files can be separated by space
rna = "050208_mRNA30minx4_RH6.cel"
## rna = "050209_mRNAx4_30min_re-hybe_RH6.cel"
## rna = "050218_polyA-RNA_RH6_4x15min.cel"

dna = "09_11_04_S96_genDNA_16hrs_45C_noDMSO.cel  041119_S96genDNA_re-hybe.cel   041120_S96genDNA_re-hybe.cel"

## Directory Path to the CEL files
celpath = "/home/huber/Marina"
#pathrna = "/ebi/research/huber/Projects/tilingArray/Celfiles"

##

  
## Now we are done with setting the parameters, and can run the code
## ------------------------------------------------------------------------
library("tilingArray")

dna = strsplit(dna, " +")[[1]]

## Step 1: read the cel files
a = readCel2eSet(c(rna, dna), path=celpath)

## Step 2: DNA normalization
library("davidTiling")  ## Get this from Bioconductor
data("probeAnno")
data("gff")

isRNA = a$filename %in% rna
isDNA = a$filename %in% dna
x = normalizeByReference(a[,isRNA], a[,isDNA], 
  pm=PMindex(probeAnno), background=BGindex(probeAnno))

## Step 3: plot along chromosome
chr=as.integer(chr)
start=as.integer(start)
end=as.integer(end)
outfilename = sprintf("%s--chr%02d:%d-%d.pdf",
   sub(".cel", "", rna, ignore.case=TRUE), chr, start, end)

pdf(file=outfilename, width=12, height=9)
plotAlongChrom(y=exprs(x), probeAnno=probeAnno, gff=gff)

dev.off()
