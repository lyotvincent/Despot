library(SAVERX)
source("sptranr/R/transpar.R")
sptFile <- "h5ads/temp1.h5spt"
sce <- Load_spt_to_SCE(sptFile)

Decontaminate_SAVERX <- function(sce){
  file <- saverx(data.matrix = sce@assays@data@listData$counts,
                 data.species = "Human", ncores = 4)
}
