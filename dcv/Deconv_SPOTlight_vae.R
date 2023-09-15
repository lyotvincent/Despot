# Title     : Using SPOTlight for Deconvolution
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/18

source("sptranr/R/transpar.R")
source("sptranr/R/_SPOTlight.R")
source("sptranr/R/_scRNA-seq.R")
library(stringr)

# decoding params
params <- fromJSON(file = "params.json")
smdFile <- params$smdFile
platform <- params$platform
if(is.null(platform)){
  platform <- "10X_Visium"  # default using 10X_Visium
}

for(decont in params$Decontamination){
  h5data <- Create_smd_h5data(decont)
  if(platform != "10X_Visium" && h5data == "SpotClean_mat"){
    message("SpotClean only Support 10X Visium data, skip it.")
    next
  }

  # read References
  sce <- Load_smdsc_to_SCE(smdFile, "sc-ref")
  rownames(sce) <- toupper(rownames(sce))
  sce <- Analysis_scRNA_seq(sce)
  sce@metadata[['HVGs']] <- toupper(rownames(sce))

  # Load the spatial transcriptome data to SE: SpatialExperiment
  spe <- Load_smd_to_SE(smdFile, h5data, platform = platform)
  # Use the symbol as rownames
  rownames(spe) <- rowData(spe)$symbol
  spotlight <- Deconvolution_SPOTlight(spe, sce)
  rownames(spotlight$mat) <- spe@colData$Barcode
  colnames(spotlight$mat) <- attr(table(sce$free_annotation), 'names')
  Save_smd_from_SPOTlight(smdFile, h5data, spotlight, pp_mtd = "vae")
}
