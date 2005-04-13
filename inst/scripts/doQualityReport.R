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
  c("041120_S96genDNA_re-hybe.cel.gz"))

hybeType=c(rep("Reverse", 3))  ## "Direct"

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

