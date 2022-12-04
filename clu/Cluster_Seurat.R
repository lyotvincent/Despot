# Title     : Using Seurat for clustering
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/10

source("sptranr/R/transpar.R")
source("sptranr/R/_Seurat.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
imgdir <- paste0(params$dataPath, "/spatial")

for(decont in params$Decontamination){
  h5data <- "matrix"
  if(decont == "SpotClean"){
    h5data <- "SpotClean_mat"
  }
  if(decont == "SPCS"){
    h5data <- "SPCS_mat"
  }

  seu <- Load_spt_to_Seurat(sptFile,imgdir, h5data)
  seu <- Preprocess_Seurat(seu)
  seu <- Cluster_Seurat(seu)
  Save_spt_from_Seurat(sptFile, seu, h5data)
}
#p1 <- DimPlot(seu, reduction = "umap", label = TRUE)
#p2 <- SpatialDimPlot(seu, label = TRUE, label.size = 3)
#pdf('spotClean_Seurat1.pdf', 10  1, 5)
#plot_grid(p1, p2)
#dev.off()
# sptFile <- "h5ads/CID4971_minor.h5spt"
# seu <- Load_sptsc_to_Seurat(sptFile)
# seu0 <- seu[, which(seu$orig.ident == "T cells CD8+")]
# seu0 <- SCTransform(seu0, verbose = F)
# resls <- c(0.5, 0.6, 0.7, 0.8, 0.9)
# seu0 <- FindVariableFeatures(seu0, nfeatures = 2000)
# for(resl in resls){
#   pdf(paste0("subcluster/CD8T_res", as.character(resl), ".pdf"), width = 5, height = 18)
#   for(dims in 11:20){
#     seu0 <- RunPCA(seu0, assay = "SCT", features = VariableFeatures(seu0))
#     seu0 <- FindNeighbors(seu0, reduction = "pca", dims = 1:dims, features = VariableFeatures(seu0))
#     seu0 <- FindClusters(seu0, verbose = FALSE, resolution = resl)
#     seu0 <- RunUMAP(seu0, reduction = "pca", dims = 1:dims)
#     markers <- FindAllMarkers(seu0, features = VariableFeatures(seu0))
#     top20 <- c()
#     for(i in 0:16){
#       marker <- markers[which(markers$cluster == i), "gene"]
#       top20 <- cbind(top20,head(marker, 15))
#     }
#     top20 <- as.character(top20)
#     plot(DoHeatmap(seu0, top20) + NoLegend())
#   }
#   dev.off()
# }
#
# pdf(paste0("subcluster/CD8T_res0.pdf"), width = 5, height = 10)
# for(dims in 5){
#   seu0 <- RunPCA(seu0, assay = "SCT")
#   seu0 <- FindNeighbors(seu0, reduction = "pca", dims = 1:dims)
#   seu0 <- FindClusters(seu0, verbose = FALSE, resolution = 0.6)
#   seu0 <- RunUMAP(seu0, reduction = "pca", dims = 1:dims)
#   markers <- FindAllMarkers(seu0)
#   top20 <- c()
#   for(i in 0:16){
#     marker <- markers[which(markers$cluster == i), "gene"]
#     top20 <- cbind(top20,head(marker, 10))
#   }
#   top20 <- as.character(top20)
#   CD8Tidents <- c("GZMB+ Tex","GZMK+ Tex", "ZNF683+ Tm", "DUSP4+ nTreg", "CD200+ terminal Tex", "TCF7+ Tn1", "iTreg", "Tn2")
#   names(CD8Tidents) <- levels(seu0)
#   seu0 <- RenameIdents(seu0, CD8Tidents)
#   plot(DoHeatmap(seu0, top20) + NoLegend())
# }
# dev.off()
# pdf(paste0("subcluster/CD8T_umap0.pdf"), width = 5, height = 5)
# plot(UMAPPlot(seu0))
#  dev.off()
#
# pdf(paste0("subcluster/CD8T_feature1.pdf"), width = 7, height = 10)
# plot(FeaturePlot(seu0, features = c("CD8A", "CD8B", "NKG7", "CD3D","CD3E","CD3G")))
# dev.off()
# sce <- as(seu0, "SingleCellExperiment")
# Save_Seurat_to_tsv("data/human_breast_cancer_atlas/GSE176078_RAW/GSM5354531_CID44971", seu0)

