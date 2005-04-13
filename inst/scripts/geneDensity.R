## calculate the gene density of different chromosomes
if(!exists("probeAnno"))load("probeAnno.rda")

sgff = gff[ which(gff$strand %in% c("-", "+")), ]

splitBy = list(sgff$seqname, sgff$strand)
ends  = split(sgff$end, splitBy)
genes = split(sgff$feature, splitBy)
  
len  = sapply(ends, max)
nrg  = sapply(genes, function(x) sum(x=="gene"))
  
gd  = nrg / len * 1e6
