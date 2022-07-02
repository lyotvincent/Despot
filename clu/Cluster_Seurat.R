# Title     : Using Seurat for clustering
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/10

source("sptranr/R/transpar.R")
source("sptranr/R/_Seurat.R")

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
  imgdir <- paste0(params$dataPath, "/spatial")

  seu <- Load_spt_to_Seurat(sptFile,imgdir, h5data)
  seu <- Preprocess_Seurat(seu)
  seu <- Cluster_Seurat(seu)
  Save_spt_from_Seurat(sptFile, seu, h5data)
}
#p1 <- DimPlot(seu, reduction = "umap", label = TRUE)
#p2 <- SpatialDimPlot(seu, label = TRUE, label.size = 3)
#pdf('spotClean_Seurat1.pdf', 10, 5)
#plot_grid(p1, p2)
#dev.off()
