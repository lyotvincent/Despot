
source("sptranr/R/transpar.R")
source("sptranr/R/_Seurat.R")
source("sptranr/R/_scRNA-seq.R")

sptFile <- "h5ads/CID4971_minor.h5spt"
seu <- Load_smdsc_to_Seurat(sptFile)
seu0 <- seu[, which(seu$orig.ident == "T cells CD8+")]
seu0 <- SCTransform(seu0, verbose = F)
resls <- c(0.5, 0.6, 0.7, 0.8, 0.9)
# seu0 <- FindVariableFeatures(seu0, nfeatures = 2000)
for(resl in resls){
  pdf(paste0("subcluster/CD8T_res", as.character(resl), ".pdf"), width = 5, height = 18)
  for(dims in 11:20){
    seu0 <- RunPCA(seu0, assay = "SCT", features = VariableFeatures(seu0))
    seu0 <- FindNeighbors(seu0, reduction = "pca", dims = 1:dims, features = VariableFeatures(seu0))
    seu0 <- FindClusters(seu0, verbose = FALSE, resolution = resl)
    seu0 <- RunUMAP(seu0, reduction = "pca", dims = 1:dims)
    markers <- FindAllMarkers(seu0, features = VariableFeatures(seu0))
    top20 <- c()
    for(i in 0:16){
      marker <- markers[which(markers$cluster == i), "gene"]
      top20 <- cbind(top20,head(marker, 15))
    }
    top20 <- as.character(top20)
    plot(DoHeatmap(seu0, top20) + NoLegend())
  }
  dev.off()
}

svg(paste0("subcluster/CD8T_res0.svg"), width = 5, height = 10)
for(dims in 6){
  seu0 <- RunPCA(seu0, assay = "SCT")
  seu0 <- FindNeighbors(seu0, reduction = "pca", dims = 1:5)
  seu0 <- FindClusters(seu0, verbose = FALSE, resolution = 0.6)
  seu0 <- RunUMAP(seu0, reduction = "pca", dims = 1:5)
  markers <- FindAllMarkers(seu0,only.pos = T)
  top20 <- c()
  for(i in 0:8){
    marker <- markers[which(markers$cluster == i), "gene"]
    top20 <- cbind(top20,head(marker, 9))
  }
  top20 <- as.character(top20)
  plot(DoHeatmap(seu0, top20) + NoLegend())
}
dev.off()
CD8Tidents <- c("GZMB+ Tex","GZMK+ Tex", "ZNF683+ Tm", "ZFP36+ Tm", "terminal Tex", "Tn1", "ENTPD1+ Tex", "Tn2")
names(CD8Tidents) <- levels(seu0)
seu0 <- RenameIdents(seu0, CD8Tidents)
svg(paste0("subcluster/CD8T_umap1.svg"), width = 6.5, height = 5)
plot(UMAPPlot(seu0))
dev.off()

# dots = c('GZMB', 'CXCL13', 'APOBEC3G', 'APOBEC3C','LAG3','NKG7',
#          'GZMH', 'CCL5','GZMK', 'GZMA', 'TRBV6-5', 'TRAV13-1',
#          'ZNF683', 'CXCR6', 'JUN', 'HOPX',
#          "TSC22D3", "ZNF331","CREM", "ZFP36","SRGN","ZFP36L2","CXCR4",
#          'CD200', 'TIGIT', 'IFNG','CCL3', 'CCL4',  'TNFRSF9',
#          'IL7R', 'TCF7', 'TC2N', 'LTB', 'RPL3', 'RPS3',
#          'TNFRSF18', 'ENTPD1', 'LAIR2', 'CD7', 'RGS13', 'HTRA1',
#          'IFIT1', 'IFIT3', 'LAMP3', 'STAT1', 'ISG15', 'MX1', 'CMPK2')

dots1 = c('CXCL13', 'GZMB',  'APOBEC3G', 'APOBEC3C',
          'GZMA','GZMK', 'TRAV13-1', 'TRBV6-5',
         'ZNF683', 'CXCR6',  'HOPX',
         "TSC22D3", "ZFP36","CXCR4", "ZNF331",
         'CCL4', 'TIGIT', 'CCL3','IFNG', 'CD200',
         'IL7R','LTB','TC2N','TCF7',
         'TNFRSF18', 'ENTPD1', 'LAIR2', 'HTRA1',
         'ISG15', 'MX1','STAT1','IFIT1',  'LAMP3')
svg(paste0("subcluster/CD8T_dot1.svg"), width = 8.5, height = 4)
plot(DotPlot(seu0, features=dots1, cols=c('lightgrey', '#800000')) +
       ggtitle("highly variable genes of CD8+ T-cell subtypes") +
       theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
             plot.title = element_text(hjust=0.5,face = 'plain')))
dev.off()

#
# pdf(paste0("subcluster/CD8T_feature1.pdf"), width = 7, height = 10)
# plot(FeaturePlot(seu0, features = c("CD8A", "CD8B", "NKG7", "CD3D","CD3E","CD3G")))
# dev.off()
# sce <- as(seu0, "SingleCellExperiment")
# Save_Seurat_to_tsv("data/human_breast_cancer_atlas/GSE176078_RAW/GSM5354531_CID44971", seu0)
sptFile <- "h5ads/V1_Adult_Mouse_Brain.h5spt"
seu <- Load_sptsc_to_Seurat(sptFile)
seu@active.ident <- as.factor(seu$orig.ident)
seu <- SCTransform(seu, verbose = F)
seu <- FindVariableFeatures(seu, nfeatures = 10000)
seu <- RunPCA(seu, assay = "SCT", features = VariableFeatures(seu))
dims <- 6
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:dims, features = VariableFeatures(seu))
seu <- RunUMAP(seu, reduction = "pca", dims = 1:dims)
UMAPPlot(seu)
markers <- FindAllMarkers(seu,only.pos = T, features = VariableFeatures(seu))
top20 <- c()
for(i in unique(seu@active.ident)){
  marker <- markers[which(markers$cluster == i), "gene"]
  top20 <- cbind(top20,head(marker, 9))
}
top20 <- as.character(top20)
plot(DoHeatmap(seu, top20) + NoLegend())
dots1 = c('AQP4', 'LCAT', 'MLC1', 'ACSBG1','HOPX',
          'ESAM','TM4SF1','FLT1','ELTD1','PECAM1',
          'NAP1L5','GAD1', 'GAD2',  'SLC32A1', 'ARX',
          'C1QA', 'C1QB', 'C1QC', 'TYROBP','PF4',
          'MBP', 'MOBP','MOG','APOD',  'PLP1',
          'HPCA', 'GRIA1', 'CPNE6', 'TSPAN13','WIPF3',
          'MEF2C', 'PCSK2', 'NECAB3', 'TBR1', 'DPP10')
svg(paste0("figures/heatmap_brain.svg"), width = 10, height = 4.5)
plot(DotPlot(seu, features=dots1, cols=c('lightgrey', '#800000')) +
       ggtitle("molecular markers of Zeisel's profiles") +
       theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
             plot.title = element_text(hjust=0.5,face = 'plain')))
dev.off()


sptFile <- "h5ads/FFPE_Mouse_Kidney_DropSeq.h5spt"
seu <- Load_sptsc_to_Seurat(sptFile)
seu@active.ident <- as.factor(seu$orig.ident)
seu <- SCTransform(seu, verbose = F)
seu <- FindVariableFeatures(seu, nfeatures = 10000)
seu <- RunPCA(seu, assay = "SCT", features = VariableFeatures(seu))
dims <- 6
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:dims, features = VariableFeatures(seu))
seu <- RunUMAP(seu, reduction = "pca", dims = 1:dims)
UMAPPlot(seu)
markers <- FindAllMarkers(seu,only.pos = T, features = VariableFeatures(seu))
top20 <- c()
for(i in unique(seu@active.ident)){
  marker <- markers[which(markers$cluster == i), "gene"]
  top20 <- cbind(top20,head(marker, 6))
}
top20 <- as.character(top20)
plot(DoHeatmap(seu, top20) + NoLegend())

dots1 = c('AQP2', 'SLC8A1', 'SLC16A7','SLC12A1', 'EPHA7', 'SLC5A12', 'SLC7A13')
# dots1 = c('AQP4', 'LCAT', 'MLC1', 'ACSBG1','HOPX',
#           'ESAM','TM4SF1','FLT1','ELTD1','PECAM1',
#           'NAP1L5','GAD1', 'GAD2',  'SLC32A1', 'ARX',
#           'C1QA', 'C1QB', 'C1QC', 'TYROBP','PF4',
#           'MBP', 'MOBP','MOG','APOD',  'PLP1',
#           'HPCA', 'GRIA1', 'CPNE6', 'TSPAN13','WIPF3',
#           'MEF2C', 'PCSK2', 'NECAB3', 'TBR1', 'DPP10')
svg(paste0("figures/heatmap_kidney_sc.svg"), width = 6, height = 4)
seu0 <- seu[, which(seu$orig.ident %in% c('CD-PC', 'CNT', 'DCT', 'LH(AL)', 'LH(DL)', 'PT(S1-S2)', 'PT(S3)'))]
plot(DotPlot(seu0, features=dots1, cols=c('lightgrey', '#800000')) +
       ggtitle("molecular markers of single-cell profiles") +
       theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
             plot.title = element_text(hjust=0.5,face = 'plain')))
dev.off()


sptFile <- "h5ads/FFPE_Mouse_Kidney_sNucSeq.h5spt"
seu <- Load_sptsc_to_Seurat(sptFile)
seu@active.ident <- as.factor(seu$orig.ident)
seu <- SCTransform(seu, verbose = F)
seu <- FindVariableFeatures(seu, nfeatures = 10000)
seu <- RunPCA(seu, assay = "SCT", features = VariableFeatures(seu))
dims <- 6
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:dims, features = VariableFeatures(seu))
seu <- RunUMAP(seu, reduction = "pca", dims = 1:dims)
UMAPPlot(seu)
markers <- FindAllMarkers(seu,only.pos = T, features = VariableFeatures(seu))
top20 <- c()
for(i in unique(seu@active.ident)){
  marker <- markers[which(markers$cluster == i), "gene"]
  top20 <- cbind(top20,head(marker, 6))
}
top20 <- as.character(top20)
plot(DoHeatmap(seu, top20) + NoLegend())

dots1 = c('AQP2','MGAT4C','FXYD4',
          'SLC8A1','SCNN1B','SCNN1G',
          'SLC12A3','TRPM6','SLC16A7',
          'EMCN','FLT1', 'ADGRL4',
          'KIT','CLNK', 'AQP6',
          'SLC26A4', 'ATP6V1C2', 'LSAMP',
          'SLC12A1', 'ERBB4','EGF',
          'NCAM1','SORCS3','EPHA7',
          'CFH','FHL2',  'AK5',
          'PTPRC','MYO1F','LST1',
          'THSD7A','NPHS1','PTPRO',
          'SLC34A1','NOX4','SLC5A12',
          'GHR','SLCO1A6','SLC7A13')
svg(paste0("figures/heatmap_kidney_sn.svg"), width = 10, height = 6)
# seu0 <- seu[, which(seu$orig.ident %in% c('CD-PC', 'CNT', 'DCT', 'LH(AL)', 'LH(DL)', 'PT(S1-S2)', 'PT(S3)'))]
plot(DotPlot(seu, features=dots1, cols=c('lightgrey', '#800000')) +
       ggtitle("molecular markers of single-nucleus profiles") +
       theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
             plot.title = element_text(hjust=0.5,face = 'plain')))
dev.off()


sptFile <- "h5ads/PDAC-A.h5spt"
seu <- Load_sptsc_to_Seurat(sptFile)
seu@active.ident <- as.factor(seu$orig.ident)
seu <- SCTransform(seu, verbose = F)
seu <- FindVariableFeatures(seu, nfeatures = 10000)
seu <- RunPCA(seu, assay = "SCT", features = VariableFeatures(seu))
dims <- 6
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:dims, features = VariableFeatures(seu))
seu <- RunUMAP(seu, reduction = "pca", dims = 1:dims)
UMAPPlot(seu)
markers <- FindAllMarkers(seu,only.pos = T, features = VariableFeatures(seu))
top20 <- c()
for(i in unique(seu@active.ident)){
  marker <- markers[which(markers$cluster == i), "gene"]
  top20 <- cbind(top20,head(marker, 10))
}
# 'Acinar cells', 'Cancer clone A', 'Cancer clone B', 'Fibroblasts'
seu0 <- seu[, which(seu$orig.ident %in% c(
  'pDCs','RBCs','T cells & NK cells','Tuft cells'
                                          ))]
top20 <- as.character(top20)
svg(paste0("figures/heatmap_PDAC.svg"), width = 10, height = 20)
plot(DoHeatmap(seu0, top20) + NoLegend())
dev.off()
dots1 = c('CELA3B','CTRB1','CTRB2',
          'TM4SF1','S100A14','CLDN1',
          'LY6D','S100A4','KRT16',
          'APOL1','SOX4','DUOX2',
          'CRISP3','AQP3','TCN1',
          'TM4SF4', 'C3', 'APCS',
          'TFF1', 'TFF2', 'TFF3',
          'STMN2','BEX1', 'CPE',
          'EMCN', 'CDH5','SOX18',
          'COL1A2', 'COL3A1', 'COL1A1',
          'FPR3','PLA2G7', 'MMP19',
          'CD68','CD52','FCN1',
          'TPSAB1', 'TPSB2', 'CPA3',
          'C1QA','C1QB', 'C1QC',
          'S100A9', 'S100A8', 'S100A12',
          'IRF4','PLD4','JCHAIN',
          'HBD', 'HBA1', 'HBB',
          'IL32','CCL5',  'NKG7',
          'GNLY','CD3D','SRGN',
          'BMX','SH2D6', 'TRPM5'
          )
svg(paste0("figures/heatmap_PDAC-A.svg"), width = 15, height = 7)
seu0 <- seu[, which(seu$orig.ident %in% c(
  'Acinar cells', 'Cancer clone A', 'Cancer clone B',
  'Ductal - APOL1 high/hypoxic',
  'Ductal - CRISP3 high/centroacinar like',
  'Ductal - MHC Class II',
  'Ductal - terminal ductal like',
  'Fibroblasts',
  'Endocrine cells',
  'Endothelial cells',
  'Macrophages A',
  'Macrophages B',
  'Mast cells',
  'mDCs A',
  'mDCs B',
  'Monocytes',
  'pDCs',
  'RBCs',
  'T cells & NK cells',
  'Tuft cells'
))]
plot(DotPlot(seu0, features=dots1, cols=c('lightgrey', '#800000')) +
       ggtitle("molecular markers of PDAC-A single-cell profiles") +
       theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
             plot.title = element_text(hjust=0.5,face = 'plain')))
dev.off()

sptFile <- "h5ads/CID4971.h5spt"
seu <- Load_sptsc_to_Seurat(sptFile)
seu@active.ident <- as.factor(seu$orig.ident)
seu <- SCTransform(seu, verbose = F)
seu <- FindVariableFeatures(seu, nfeatures = 10000)
seu <- RunPCA(seu, assay = "SCT", features = VariableFeatures(seu))
dims <- 6
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:dims, features = VariableFeatures(seu))
seu <- RunUMAP(seu, reduction = "pca", dims = 1:dims)
UMAPPlot(seu)
markers <- FindAllMarkers(seu,only.pos = T, features = VariableFeatures(seu))
top20 <- c()
for(i in unique(seu@active.ident)){
  marker <- markers[which(markers$cluster == i), "gene"]
  top20 <- cbind(top20,head(marker, 10))
}
top20 <- as.character(top20)
plot(DoHeatmap(seu, top20) + NoLegend())
dots1 = c("MS4A1","CD79A","NCF1","CD79B","CD19",
          "APOD", "DCN", "LUM","COL1A1","CTGF",
          "MGST1","MARCKSL1", "CD24", "UCHL1","STMN1",
          "RNASE1","PLVAP",'PECAM1',"RBP7" ,"FABP4",
          'CST3',"TYROBP",'C1QA','C1QB','C1QC',
          "AZGP1","CXCL2","LTF","EPCAM","SCGB3A1",
          "JCHAIN","IGLL5","IGHG4","IGHG3","IGHG2",
          "NDUFA4L2","RGS5","ACTA2","GJA4","HIGD1B",
          "IL32", "CD3E","CD7","CCL5","NKG7"
)
svg(paste0("figures/heatmap_CID4971.svg"), width = 10.5, height = 5)
# seu0 <- seu[, which(seu$orig.ident %in% c('CD-PC', 'CNT', 'DCT', 'LH(AL)', 'LH(DL)', 'PT(S1-S2)', 'PT(S3)'))]
plot(DotPlot(seu, features=dots1, cols=c('lightgrey', '#800000')) +
       ggtitle("molecular markers of CID4971 single-cell profiles") +
       theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
             plot.title = element_text(hjust=0.5,face = 'plain')))
dev.off()
