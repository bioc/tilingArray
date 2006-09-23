## A script for Marina Granovskaia to make Along-Chromosome Plots
## (C) W. Huber 2006
##
## Please specify the parameters in the following lines:

## Name of RNA and of DNA hybes. Multiple files can be separated by comma.
rna = "050208_mRNA30minx4_RH6.cel, 050209_mRNAx4_30min_re-hybe_RH6.cel, 050218_polyA-RNA_RH6_4x15'.cel"
dna = "09_11_04_S96_genDNA_16hrs_45C_noDMSO.cel, 041119_S96genDNA_re-hybe.cel, 041120_S96genDNA_re-hybe.cel"

## Directory Path to the CEL files
celpath = "/home/huber/Marina/Celfiles"
#pathrna = "/ebi/research/huber/Projects/tilingArray/Celfiles"

  
## Now we are done with setting the parameters, and can run the code
## ------------------------------------------------------------------------

## Step 1: read the cel files "a"
a <- readCel2eSet(c(rna, dna), path=celpath)
