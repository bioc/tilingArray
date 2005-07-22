
library("tilingArray")
library("affy")

options(error=recover)
source("~/madman/Rpacks/tilingArray/R/qualityReport.R")

if(!exists("probeAnno"))
  load("probeAnno.rda")

files = file.path("Celfiles",
  c("050621_dirPolyARNA_10ug_2-3.cel.gz",
    "050621_dirPolyARNA_10ug_2-3_4x.cel.gz"))

normRefFiles = file.path("Celfiles",
  c("09_11_04_S96_genDNA_16hrs_45C_noDMSO.cel.gz",
    "041119_S96genDNA_re-hybe.cel.gz",
    "041120_S96genDNA_re-hybe.cel.gz"))

hybeType=c("Reverse", "Direct")[2]

if(!exists("x"))
  x = read.affybatch(filenames=files, compress=TRUE, verbose=TRUE)

if(!exists("normRef"))
  normRef = read.affybatch(filenames=normRefFiles, compress=TRUE, verbose=TRUE)

for(i in seq(along=files))
  qualityReport(x = x[,i], hybeType=hybeType, normRef = normRef,
              selectGenes = c("YLR342W", "YOR153W", "YDL185W", "YKL182W",
                "YGR027C", "YDR384C"),
              gff = gff, probeAnno = probeAnno,
              outputDir = "qualityReports")

