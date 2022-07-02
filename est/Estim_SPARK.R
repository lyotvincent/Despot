# Title     : Using SPARK for Estimation high variable genes
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/12

source("sptranr/R/transpar.R")
source("sptranr/R/_SPARK.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
for(decont in params$Decontamination){
  h5data <- "matrix"
  if(decont == "SpotClean"){
    h5data <- "SpotClean_mat"
  }
  if(decont == "SPCS"){
    h5data <- "SPCS_mat"
  }

  sce <- Load_spt_to_SCE(sptFile, h5data = h5data)
  if(h5data == 'matrix'){
    rownames(sce) <- sce@rowRanges@elementMetadata@listData$gene_id
  }else{
    rownames(sce) <- sce@rowRanges@elementMetadata@listData$gene_name
  }
  sparkX <- Estimation_SPARKX(sce)

  Save_spt_from_SPARK(sptFile, sparkX, h5data)

  print(paste0("Estimation with `SPARK` finished, HVGs saved in /", h5data, '/features/is_HVG/SPARK'))
}
