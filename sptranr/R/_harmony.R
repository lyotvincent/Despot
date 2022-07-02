library(harmony)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(celldex)
library(scater)
library(SingleR)
library(BiocParallel)
hpca.se <- celldex::HumanPrimaryCellAtlasData()


root <- "/home/rzh/PycharmProjects/spatial/data/scRNA-seq/breast cancer/nju"

dts <- list.files(root)

scQC <- function(seu,lg=100,hg=100000,lm=-Inf,hm=0.1){
  #QC and filter Seurat object

  #Inputs:
  #seu = Seruat object

  #Calculate %mito genes
  mito.genes <- grep(pattern="^MT-",x=rownames(seu),value=TRUE)
  if(length(mito.genes > 0)){
    percent.mito <- colSums(seu[mito.genes,])/colSums(seu)
    seu <- AddMetaData(object=seu,metadata=percent.mito,col.name="percent.mito")


    counts.cell.raw <- dim(seu)[2]
    counts.gene.raw <- dim(seu)[1]

    seu <- seu[, seu[["percent.mito"]] < hm & lg < seu[['nFeature_RNA']] & seu[['nFeature_RNA']] < hg]
  }

  #Plot filtered stats
  counts.cell.filtered <- dim(seu)[2]

  #Remove mito genes
  seu <- seu[!rownames(seu) %in% mito.genes,]

  counts.gene.filtered <- dim(seu)[1]

  return(seu)
}

exp1 <- read.table(paste0(datapath, '/', dts[4]), header = T, sep = '\t')
exp1 <- Matrix(data = as.matrix(exp1))
exp2 <- read.table(paste0(datapath, '/', dts[5]), header = T, sep = '\t')
exp2 <- Matrix(data = as.matrix(exp2))
exp3 <- read.table(paste0(datapath, '/', dts[6]), header = T, sep = '\t')
exp3 <- Matrix(data = as.matrix(exp3))

genes <- Reduce(intersect ,list(rownames(exp1), rownames(exp2), rownames(exp3)))
exp1 <- exp1[genes,]
exp2 <- exp2[genes,]
exp3 <- exp3[genes,]

seu1 <- CreateSeuratObject(exp1)
Idents(seu1) <- "P2A"
seu2 <- CreateSeuratObject(exp2)
Idents(seu2) <- "P2B"
seu3 <- CreateSeuratObject(exp3)
Idents(seu3) <- "P2C"

seu <- Reduce(function(x,y) merge(x,y,all=T),c(seu1, seu2, seu3))
seu$orig.ident <- seu@active.ident
seu <- scQC(seu)
seu <- Seurat::NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

seu <- seu %>% RunHarmony("orig.ident", plot_convergence = TRUE)
seu <- seu %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.5) %>%
  identity()
seu <- seu %>%
  RunTSNE(reduction = "harmony", dims = 1:30)

logFCfilter=0.5
adjPvalFilter=0.05
seu.markers <- FindAllMarkers(object = seu,
                                only.pos = FALSE,
                                min.pct = 0.25,
                                logfc.threshold = logFCfilter)
sig.markers=seu.markers[(abs(as.numeric(as.vector(seu.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(seu.markers$p_val_adj))<adjPvalFilter),]

common_hpca <- intersect(rownames(seu), rownames(hpca.se))
hpca.se <- hpca.se[common_hpca,]
red.seu <- SingleR(test = seu@assays$RNA@data,
                        ref = hpca.se,
                        labels = hpca.se$label.main,
                        clusters = seu@active.ident,
                        fine.tune = TRUE, BPPARAM = MulticoreParam(40))

celltype <- as.numeric(Idents(seu))
lvs <- unique(celltype)
for(lv in lvs){
  ct <- red.seu$labels[lv]
  celltype[which(celltype == lv)] <- ct
}
celltype <- as.factor(celltype)
seu@meta.data$free_annotation <- celltype
