source("sptranr/R/transpar.R")
source("sptranr/R/_Seurat.R")
source("sptranr/R/_scRNA-seq.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile

seu.sc <- Load_sptsc_to_Seurat(sptFile, h5data = "scRNA_seq")
seu.sc <- Preprocess_Seurat(seu.sc, assay = "RNA")

for(decont in params$Decontamination){
  h5data <- "matrix"
  if(decont == "SpotClean"){
    h5data <- "SpotClean_mat"
  }
  if(decont == "SPCS"){
    h5data <- "SPCS_mat"
  }
  imgdir <- paste0(params$dataPath, "/spatial")

  seu.sp <- Load_spt_to_Seurat(sptFile,imgdir, h5data)
  seu.sp <- Preprocess_Seurat(seu.sp, assay = "Spatial")
  predictions <- Deconvolution_Seurat(seu.sp, seu.sc)
  Save_spt_from_Seurat.dcv(sptFile, h5data, predictions)
}
