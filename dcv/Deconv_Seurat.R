source("sptranr/R/transpar.R")
source("sptranr/R/_Seurat.R")
source("sptranr/R/_scRNA-seq.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
params <- h5read(sptFile, "configs")
imgdir <- paste0(params$dataPath, "/spatial")
platform <- params$platform
if(is.null(platform)){
  platform <- "10X_Visium"  # default using 10X_Visium
}

seu.sc <- Load_sptsc_to_Seurat(sptFile, h5data = "scRNA_seq")
seu.sc <- Preprocess_Seurat(seu.sc, assay = "RNA")

for(decont in params$Decontamination){
  h5data <- Create_spt_h5data(decont)
  if(platform != "10X_Visium" && h5data == "SpotClean_mat"){
    message("SpotClean only Support 10X Visium data, skip it.")
    next
  }

  seu.sp <- Load_spt_to_Seurat(sptFile,
                               platform = platform,
                               imgdir = imgdir,
                               h5data = h5data)
  seu.sp <- Preprocess_Seurat(seu.sp, assay = "Spatial")
  predictions <- Deconvolution_Seurat(seu.sp, seu.sc)
  Save_spt_from_Seurat.dcv(sptFile, h5data, predictions)
}
