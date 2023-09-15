library(SAVERX)
source("sptranr/R/transpar.R")
smdFile <- "h5ads/temp1.h5spt"
sce <- Load_smd_to_SCE(smdFile)

Decontaminate_SAVERX <- function(sce){
  file <- saverx(data.matrix = sce@assays@data@listData$counts,
                 data.species = "Human", ncores = 4)
}
