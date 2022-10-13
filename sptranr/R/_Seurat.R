library(Seurat)
library(SeuratObject)
library(cowplot)
library(Matrix)
library(Rcpp)

Preprocess_Seurat <- function(seu, assay = "Spatial"){
  seu <- SCTransform(seu, assay = assay,
                         return.only.var.genes = FALSE, verbose = FALSE)
  seu <- RunPCA(seu, assay = "SCT", verbose = FALSE)
  return(seu)
}

Cluster_Seurat <- function(seu){
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30)
  seu <- FindClusters(seu, verbose = FALSE)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:30)
  return(seu)
}

Estimation_Seurat <- function(seu){
   de_markers <-FindAllMarkers(seu, assay = "SCT")
   return(de_markers)
}

Deconvolution_Seurat <- function(seu.sp, seu.sc){
  anchors <- FindTransferAnchors(reference = seu.sc, query = seu.sp, normalization.method = "SCT")
  predictions.assay <- TransferData(anchorset = anchors, refdata = seu.sc$orig.ident, prediction.assay = TRUE,
                                    weight.reduction = seu.sp[["pca"]], dims = 1:30)
  return(predictions.assay)
}

FindMarkers_Seurat <- function(seu, ident){
  de_markers <- FindMarkers(seu, ident.1 = ident)
  return(de_markers)
}

Save_spt_from_Seurat.dcv <- function(sptFile, h5data, predictions.assay){

  results <- Matrix(t(predictions.assay@data[-dim(predictions.assay@data)[1],]))
  results <- as(results, "dgeMatrix")
  cell_type_names <- rownames(predictions.assay)
  h5createGroup(sptFile, paste0(h5data, '/deconv/Seurat'))

  # save 1d weights
  Create_spt_array1d(sptFile,
                     arr = results@x,
                     sptloc = paste0(h5data, '/deconv/Seurat/weights'),
                     mode = 'double')
  # save shape
  Create_spt_array1d(sptFile,
                     arr = results@Dim,
                     sptloc = paste0(h5data, '/deconv/Seurat/shape'),
                     mode = 'integer')
  # save dim names
  Create_spt_array1d(sptFile,
                     arr = rownames(results),
                     sptloc = paste0(h5data, '/deconv/Seurat/barcodes'),
                     mode = 'character')
  Create_spt_array1d(sptFile,
                     arr = colnames(results),
                     sptloc = paste0(h5data, '/deconv/Seurat/cell_type'),
                     mode = 'character')
}

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

  #Filter/normalize data
  seu <- NormalizeData(object=seu, normalization.method="LogNormalize",scale.factor=10000)

  #Plot filtered stats
  counts.cell.filtered <- dim(seu)[2]

  #Remove mito genes
  seu <- seu[!rownames(seu) %in% mito.genes,]

  counts.gene.filtered <- dim(seu)[1]

  return(seu)
}

Load_txtsc_to_Seurat <- function(txt.list, root = "", header = T, sep = '\t'){
  data.list <- list()
  gene.list <- list()
  # load datas from txt.list to exps
  for(mat in txt.list){
    datapath <- paste0(root, '/', mat)
    message(paste0("Load scRNA-seq from ", datapath))
    exp <- read.table(datapath, header = header, sep = sep)
    # store with sparse dgCMatrix
    data.list[[mat]] <- Matrix(data = as.matrix(exp))
    gene.list[[mat]] <- rownames(data.list[[mat]])
  }
  rm(exp)

  # get intersection among different datas
  genes <- Reduce(intersect ,gene.list)
  rm(gene.list)

  message("Merging scRNA-seq data...")
  seu.list <- list()
  library(SeuratObject)
  # remove genes out of intersection and Create Seurat Object
  for(mat in txt.list){
    seu <- CreateSeuratObject(data.list[[mat]][genes, ])
    seu@meta.data$orig.ident <- mat
    seu.list[[mat]] <- seu
  }
  rm(data.list)
  seu <- Reduce(merge, seu.list)
  rm(seu.list)
  return(seu)
}

Analysis_scRNA_seq_Seurat <- function(seu, itg.method = "harmony",
                                      anno.method = "SingleR",
                                      top.HVGs = 3000,
                                      resolution = 0.8
                                      ){
  seu <- scQC(seu)
  seu <- Seurat::NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
  seu <- ScaleData(seu)
  seu <- Seurat::RunPCA(seu, npcs = 50, verbose = FALSE)
  pc.use <- (seu@reductions$pca@stdev)^2
  pc.use <- pc.use/sum(pc.use)
  pc.use <- cumsum(pc.use)
  pc.use <- min(which(pc.use>=0.9))
  reduction <- "pca"
  if(itg.method == "harmony"){
    library(harmony)
    seu <- seu %>% RunHarmony("orig.ident")
    reduction <- "harmony"
  }
  seu <- seu %>%
    RunUMAP(reduction = reduction, dims = 1:pc.use) %>%
    FindNeighbors(reduction = reduction, dims = 1:pc.use) %>%
    FindClusters(resolution = resolution) %>%
    identity()

  logFCfilter=0.5
  adjPvalFilter=0.05
  seu.markers <- FindAllMarkers(object = seu,
                                only.pos = FALSE,
                                min.pct = 0.25,
                                logfc.threshold = logFCfilter)
  sig.markers=seu.markers[(abs(as.numeric(as.vector(seu.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(seu.markers$p_val_adj))<adjPvalFilter),]

  if(anno.method == "SingleR"){
    library(celldex)
    library(SingleR)
    library(BiocParallel)
    hpca.se <- celldex::HumanPrimaryCellAtlasData()
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
  }
  return(seu)
}

Load_10X_Visium_to_Seurat <- function(dir){
  ref <- read.csv("/home/rzh/Downloads/Breast cancer.csv", header = T)
  ref <- select(ref, "celltype", "gene")
  celltype <- unique(ref$celltype)
  markers <- list()
  for(ct in celltype){
    gene <- ref[which(ref$celltype == ct),]$gene
    markers[[ct]] <- gene
  }
  dat <- read.table(dir, sep = '\t', header = T)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat <- Seurat::Read10X(dir)
  seu <- CreateSeuratObject(dat)
  seu <- scQC(seu, lg = 20, hg = 100000)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 5000)
  seu <- ScaleData(seu)
  seu <- RunPCA(object=seu,pc.genes=VariableFeatures(seu),pcs.compute=50)
  #Calculate PCs to use
  pc.use <- (seu@reductions$pca@stdev)^2
  pc.use <- pc.use/sum(pc.use)
  pc.use <- cumsum(pc.use)
  pc.use <- min(which(pc.use>=0.9))
  seu <- RunTSNE(object=seu,reduction.use="pca",dims.use=1:pc.use)
  seu <- FindNeighbors(seu , dims = 1:pc.use)
  seu <- RunUMAP(object=seu,reduction="pca", dims = 1:pc.use)
  seu <- FindClusters(object=seu,resolution=0.2)
  seu <- BuildClusterTree(seu,reorder.numeric=TRUE)
  ##寻找差异表达的特征
  logFCfilter=0.5
  adjPvalFilter=0.05
  seu.markers <- FindAllMarkers(object = seu,
                                 only.pos = FALSE,
                                 min.pct = 0.25,
                                 logfc.threshold = logFCfilter)
  top30 <- seu.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
  top30["celltype"] <- apply(top30, 1, function(x){
                                has <- x["gene"] %in% ref$gene
                                if (has){
                                  return(ref[which(ref$gene == x["gene"]),1])
                                }else{
                                  return("unknown")
                                }})
  library(SCINA)
  library(preprocessCore)
  marker_path <- "data/scRNA-seq/prostate cancer/marker references/"
  markers <- read.csv(paste0(marker_path, "prostate markers.csv"), header = T)
  # marker_file <- list.files("data/scRNA-seq/prostate cancer/marker references")
  # markers <- list()
  # for(fn in marker_file){
  #   cell_name <- strsplit(fn, '\\.')[[1]][3]
  #   mg <- read.csv(paste0(marker_path, fn))
  #   markers[[cell_name]] <- mg$X
  # }
  exp_raw <- seu[VariableFeatures(seu),]
  exp_raw <- as.matrix(exp_raw@assays$RNA@data)
  # exp_raw=log(exp_raw+1)
  # exp_raw[]=normalize.quantiles(exp_raw)
  # exp_raw <- as.numeric(exp_raw)
  # my_markers <- in_markers[1:6]
  markers <- SCINA::preprocess.signatures(paste0(marker_path, "prostate markers.csv"))
  library(stringr)
  markers <- lapply(markers, str_to_upper)
  in_markers <- list()
  for(m in attr(markers, "names")){
    in_marker <- markers[[m]]
    in_marker <- in_marker[in_marker %in% VariableFeatures(seu)]
    in_markers[[m]] <- in_marker
  }
  results <- SCINA::SCINA(exp_raw, markers, max_iter = 200, convergence_n = 10,
                          convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE)
  Idents(seu) <- results
  return(seu)
}

