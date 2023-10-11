source("sptranr/R/transpar.R")
source("sptranr/R/_Seurat.R")

# decoding params
params <- fromJSON(file = "params.json")
smdFile <- params$smdFile
imgdir <- paste0(params$dataPath, "/spatial")
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
  seu <- Load_smd_to_Seurat(smdFile,
                            platform = platform,
                            imgdir = imgdir,
                            h5data = h5data)
  seu <- Preprocess_Seurat(seu)
  seu <- Cluster_Seurat(seu)
  Save_smd_from_Seurat(smdFile, seu, h5data)
}


