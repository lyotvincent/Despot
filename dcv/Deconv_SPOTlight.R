# Title     : Using SPOTlight for Deconvolution
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/18

source("sptranr/R/transpar.R")
source("sptranr/R/_SPOTlight.R")
source("sptranr/R/_scRNA-seq.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
refdir <- paste0(params$tempdir, "/references")


for(decont in params$Decontamination){
  h5data <- "matrix"
  if(decont == "SpotClean"){
    h5data <- "SpotClean_mat"
  }
  if(decont == "SPCS"){
    h5data <- "SPCS_mat"
  }

  # read References
  sce <- Load_sptsc_to_SCE(sptFile, "scRNA_seq")
  rownames(sce) <- toupper(rownames(sce))
  sce <- Analysis_scRNA_seq(sce)
  # sce@metadata[['HVGs']] <- rownames(sce)

  # Load the spatial transcriptome data to SE: SpatialExperiment
  spe <- Load_spt_to_SE(sptFile, h5data)
  # Use the symbol as rownames
  rownames(spe) <- rowData(spe)$symbol
  spotlight <- Deconvolution_SPOTlight(spe, sce)
  rownames(spotlight$mat) <- spe@colData$Barcode
  colnames(spotlight$mat) <- attr(table(sce$free_annotation), 'names')
  Save_spt_from_SPOTlight(sptFile, h5data, spotlight, pp_mtd = "raw")
}
