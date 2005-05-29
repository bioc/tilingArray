
library("tilingArray")
library("affy")

options(error=recover)
source("~/madman/Rpacks/tilingArray/R/plotAlongChrom.R")
source("~/madman/Rpacks/tilingArray/R/qualityReport.R")

if(!exists("probeAnno"))
  load("probeAnno.rda")

normRefFiles = file.path("Celfiles",
  c("09_11_04_S96_genDNA_16hrs_45C_noDMSO.cel.gz",
    "041119_S96genDNA_re-hybe.cel.gz",
    "041120_S96genDNA_re-hybe.cel.gz"))

files = file.path("Celfiles",
    "050507_dirRNA_10ug_F1.cel.gz")

hybeType=c("Reverse", "Direct")[2]

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

