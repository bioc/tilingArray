library("tilingArray")
library("affy")

source("/home/huber/madman/Rpacks/tilingArray/R/qualityReport.R")
source("/home/huber/madman/Rpacks/tilingArray/R/plotAlongChrom.R")
options(error=recover)

if(!exists("probeAnno"))
  load("probeAnno.rda")

normRefFiles = file.path("Celfiles",
  c("09_11_04_S96_genDNA_16hrs_45C_noDMSO.cel.gz",
    "041119_S96genDNA_re-hybe.cel.gz",
    "041120_S96genDNA_re-hybe.cel.gz"))

files = file.path("Celfiles",
  c("050209_mRNAx4_30min_re-hybe_RH6.cel.gz",
    "041203_S96_polyAx1_RH6.cel.gz",
    "050218_polyA-RNA_RH6_4x15min.cel.gz",
    "05242_totRNA_15ugS96_dir#3.cel.gz",
    "030505_totcDNA_15ug_affy.cel.gz"))[4]

hybeType=c("Reverse", "Reverse", "Reverse", "Direct", "Reverse")[4]

if(!exists("x"))
  x = read.affybatch(filenames=files, compress=TRUE, verbose=TRUE)

if(!exists("normRef"))
  normRef = read.affybatch(filenames=normRefFiles, compress=TRUE, verbose=TRUE)

for(i in seq(along=files))
  qualityReport(x = x[,i], hybeType=hybeType[i], normRef = normRef,
              selectGenes = c("YLR342W", "YOR153W", "YDL185W", "YKL182W",
                "YGR027C", "YDR384C"),
              gff = gff, probeAnno = probeAnno,
              outputDir = "qualityReports")

