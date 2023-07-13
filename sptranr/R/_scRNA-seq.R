library(scater)
library(scran)
library(SingleCellExperiment)
library(scuttle)


Analysis_TabulaData <- function(sce, age = "18m",
                                lower.bound = 25 ,upper.bound = 2000,
                                filtered.celltype = c("nan", "CD45"),
                                top.HVGs = 3000){
  # filtering age
  sce <- sce[,sce$age %in% age]

  # filtering celltypes by bounds and the types given
  low <- table(sce$free_annotation) < lower.bound
  filtered.celltype <- append(filtered.celltype,
                              attr(low[low==T], "names"))
  up <- table(sce$free_annotation) > upper.bound
  filtered.celltype <- append(filtered.celltype,
                              attr(up[up==T], "names"))
  sce <- sce[,! sce$free_annotation %in% filtered.celltype]

  # Get the top highly variable genes
  genes <- !grepl("^Rp[l|s]|Mt", rownames(sce))
  sce <- sce[genes,]
  sce <- Analysis_scRNA_seq(sce, top.HVGs = top.HVGs)
  return(sce)
}

Analysis_scRNA_seq <- function(sce, top.HVGs = 3000){
  sce <- logNormCounts(sce)
  dec <- modelGeneVar(sce)
  hvg <- getTopHVGs(dec, n = top.HVGs)
  colLabels(sce) <- colData(sce)$free_annotation
  # store highly variable genes in sce
  sce@metadata[["HVGs"]] <- hvg

  # Get the marker genes
  mgs <- scoreMarkers(sce)
  mgs_ls <- lapply(names(mgs), function(i){
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    AUC.thresh <- 0.8
    while (AUC.thresh > 0){
      if (sum(x$mean.AUC > AUC.thresh) > 0){
        x <- x[x$mean.AUC > AUC.thresh, ]
        break
      }else{
        AUC.thresh <- AUC.thresh - 0.1
      }
    }

    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
  })
  mgs_df <- do.call(rbind, mgs_ls)
  # store marker genes in metadata
  sce@metadata[["marker_genes"]] <- mgs_df
  return(sce)
}


Analysis_scRNA_seq_noanno <- function(sce, top.HVGs = 3000,
                                      lower.bound = 25,
                                      filtered.celltype = c("nan", "unknown"),
                                      itg.method = NA,
                                      anno.signature = "data/scRNA-seq/prostate cancer/marker references/prostate cancer.csv",
                                      anno.method = "SCINA"){
  sce <- addPerCellQCMetrics(sce)
  sce <- addPerFeatureQCMetrics(sce)
  # 选取总counts数大于100的
  keep.total <- sce$total > 500
  # 选取表达的基因数大于20的
  keep.n <- sce$detected > 20
  # 根据设定的条件进行过滤
  sce <- sce[,keep.total & keep.n]

  keep_feature <- nexprs(sce, byrow=TRUE) >= 10
  sce <- sce[keep_feature,]

  sce <- logNormCounts(sce)
  dec <- modelGeneVar(sce)
  hvg <- getTopHVGs(dec, n = top.HVGs)
  sce@metadata[["HVGs"]] <- hvg
  sce <- runPCA(sce)

  reduction <- "PCA"
  if(itg.method == "harmony"){
    library(harmony)
    sce <- RunHarmony(sce, "orig.ident", plot_convergence = TRUE)
  }
  sce <- runUMAP(sce)
  # In this case, using the PCs that we chose from getClusteredPCs().
  g <- buildSNNGraph(sce, use.dimred="PCA")
  cluster <- igraph::cluster_walktrap(g)$membership

  # Usinf SCINA to generate cell annotations
  library(SCINA)
  library(preprocessCore)
  # load markers from csv.file
  markers <- SCINA::preprocess.signatures(anno.signature)
  # remove the genes not in HVGs
  in_markers <- list()
  for(m in attr(markers, "names")){
    in_marker <- markers[[m]]
    in_marker <- in_marker[in_marker %in% sce@rowRanges@elementMetadata@listData$gene_name]
    in_markers[[m]] <- in_marker
  }
  exp_raw <- as.matrix(sce@assays@data$logcounts)
  rownames(exp_raw) <- sce@rowRanges@elementMetadata@listData$gene_name
  # do SCINA
  results <- SCINA::SCINA(exp_raw, in_markers, max_iter = 200, convergence_n = 10,
                          convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE)

  # Assigning to the 'colLabels' of the 'sce'.
  colLabels(sce) <- results$cell_labels

  # filter out the Celltype with cells less than lower.bound
  low <- table(colLabels(sce)) < lower.bound
  filtered.celltype <- append(filtered.celltype,
                              attr(low[low==T], "names"))

  # filter out the cell tpye in filtered.celltype
  sce <- sce[,! colLabels(sce) %in% filtered.celltype]

  # find marker genes
  mgs <- scoreMarkers(sce)
  mgs_ls <- lapply(names(mgs), function(i){
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.8, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
  })
  mgs_df <- do.call(rbind, mgs_ls)
  # store marker genes in metadata
  sce@metadata[["marker_genes"]] <- mgs_df
  return(sce)
}


# Since the TabulaData is too large, we store the preprocessed scRNA-seq Data
Save_scRNAseq_to_spt <- function(sptFile, sce){
  # create group for scRNA-seq
  rhdf5::h5createGroup(sptFile, "scRNA_seq")

  # create DataSets for matrix
  mat <- sce@assays@data@listData$counts
  mat <- Matrix(mat, nrow = dim(mat)[1], ncol = dim(mat)[2], dimnames = dimnames(mat))

  # save 1d weights
  Create_spt_array1d(sptFile,
                     arr = mat@x,
                     sptloc = "scRNA_seq/data",
                     mode = 'integer')
  Create_spt_array1d(sptFile,
                     arr = mat@i,
                     sptloc = "scRNA_seq/indices",
                     mode = 'integer')
  Create_spt_array1d(sptFile,
                     arr = mat@p,
                     sptloc = "scRNA_seq/indptr",
                     mode = 'integer')
  Create_spt_array1d(sptFile,
                     arr = mat@Dim,
                     sptloc = "scRNA_seq/shape",
                     mode = 'integer')
  Create_spt_array1d(sptFile,
                     arr = colnames(mat),
                     sptloc = "scRNA_seq/barcodes",
                     mode = 'character')

  # create features
  h5createGroup(sptFile, "scRNA_seq/features")

  # save original names
  Create_spt_array1d(sptFile,
                     arr = toupper(rownames(mat)),
                     sptloc = "scRNA_seq/features/name",
                     mode = 'character')

  # save HVGs
  Create_spt_array1d(sptFile,
                     arr = toupper(sce@metadata$HVGs),
                     sptloc = "scRNA_seq/features/HVGs",
                     mode = 'character')

  # create idents
  h5createGroup(sptFile, "scRNA_seq/idents")

  # save annotation and index
  Create_spt_array1d(sptFile,
                     arr = sce@colData@rownames,
                     sptloc = "scRNA_seq/idents/index",
                     mode = 'character')
  Create_spt_array1d(sptFile,
                     arr = as.character(colLabels(sce)),
                     sptloc = "scRNA_seq/idents/annotation",
                     mode = 'character')
}

# some sofrwares only support .tsv files, so we provide the interface
Save_TabulaData_to_tsv <- function(outdir, sce){
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }

  # save cnt.tsv
  cnt <- sce@assays@data@listData$counts
  write.table(cnt, file = paste0(outdir, '/cnt.tsv'), sep = '\t', quote = F)
  # save mta.tsv
  mta <- data.frame(bio_celltype = as.character(sce@colData$free_annotation),
                    row.names = sce@colData@rownames)
  write.table(mta, file = paste0(outdir, '/mta.tsv'), sep = '\t', quote = F)
}

Save_Seurat_to_tsv <- function(outdir, seu0){
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  # save mtx
  cnt <- seu0@assays$SCT@counts
  writeMM(cnt, file = paste0(outdir, '/matrix.mtx'))
  # save barcodes
  bar <- colnames(seu0)
  write.table(bar, file = paste0(outdir, '/barcodes.tsv'), sep = '\t', quote = F, row.names = F, col.names = F)
  # save features
  fea <- rownames(seu0)
  write.table(fea, file = paste0(outdir, '/features.tsv'), sep = '\t', quote = F, row.names = F, col.names = F)
  # save annotation
  mta <- data.frame(orig.ident = as.character(seu0@meta.data$orig.ident),
                    bio_celltype = as.character(seu0@active.ident),
                    row.names = bar)
  write.table(mta, file = paste0(outdir, '/mta.csv'), sep = ',', quote = F)
}


Load_ref_to_SCE <- function(refdir){
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  cnt <- read.table(paste0(refdir, "/cnt_data.tsv"), header = T, sep = '\t')
  rownames(cnt) <- cnt[,1]
  cnt <- t(Matrix(as.matrix(cnt[,-1])))
  mta <- read.table(paste0(refdir, "/mta_data.tsv"), header = T, sep = '\t')
  mta <- as.factor(mta[,-1])
  se <- SummarizedExperiment(assays = list(counts = cnt),
                             rowData = DataFrame(gene_name = rownames(cnt)),
                             colData = DataFrame(free_annotation = mta))
  sce <- as(se, "SingleCellExperiment")
  return(sce)
}

Load_sptsc_to_SCE <- function(sptFile, h5data = "scRNA_seq"){
  h5_obj <- rhdf5::h5read(sptFile, h5data)
  dat <- sparseMatrix(i = h5_obj$indices[] + 1,
                      p = h5_obj$indptr[],
                      x = as.numeric(h5_obj$data[]),
                      dims = h5_obj$shape[],
                      dimnames = list(NULL, NULL),
                      repr = "C")
  colData <- DataFrame(barcodes = h5_obj$barcodes,
                       free_annotation = h5_obj$idents$annotation)
  se <- SummarizedExperiment(assays = list(counts = dat),
                             rowData = DataFrame(gene_name = toupper(h5_obj$features$name)),
                             colData = colData)
  sce <- as(se, "SingleCellExperiment")
  colnames(sce) <- h5_obj$barcodes
  rownames(sce) <- toupper(h5_obj$features$name)
  return(sce)
}

Load_sptsc_to_Seurat <- function(sptFile, h5data = "scRNA_seq"){
  library(SeuratObject)
  h5_obj <- rhdf5::h5read(sptFile, h5data)
  dat <- sparseMatrix(i = h5_obj$indices[] + 1,
                      p = h5_obj$indptr[],
                      x = as.numeric(h5_obj$data[]),
                      dims = h5_obj$shape[],
                      dimnames = list(toupper(h5_obj$features$name), h5_obj$barcodes),
                      repr = "C")
  colData <- DataFrame(barcodes = h5_obj$barcodes,
                       free_annotation = h5_obj$idents$annotation)
  seu <- CreateSeuratObject(counts = dat)
  seu@meta.data$orig.ident <- h5_obj$idents$annotation
  return(seu)
}

# Load Tabula scRNA-seq data to SingleCellExperiment
Load_TabulaData_to_SCE <- function(TabulaType, TabulaData){
  library(TabulaMurisSenisData)
  listData <- listTabulaMurisSenisTissues(TabulaType)
  if(TabulaData %in% listData){
    if(TabulaType == "Droplet"){
      sce <- TabulaMurisSenisDroplet(tissues = TabulaData)[[TabulaData]]
    }else if(TabulaType == "FACS"){
      sce <- TabulaMurisSenisFACS(tissues = TabulaData)[[TabulaData]]
    }else{
      message("Tabula Type Error.")
      return(FALSE)
    }
  }else{
    message(paste0("No TabulaData:", TabulaData, " exist in ", TabulaType, "."))
    return(FALSE)
  }
  return(sce)
}

# Load 10X scRNA-seq data to SingleCellExperiment
Load_10Xsc_to_SCE <- function(tenXdir){
  library(Seurat)
  library(stringr)
  sce <- Read10X(tenXdir)
  gene_name <- rownames(sce)
  gene_name <- str_to_upper(gene_name)
  sce <- SingleCellExperiment(assays = list(counts = sce),
                              rowData = DataFrame(gene_name = gene_name),
                              colData = DataFrame(barcodes = colnames(sce)))
  return(sce)
}

# Load scRNA-seq data by the users' definition: mtx barcodes and features
Load_sc_to_SCE <- function(scbar, scfea, scmat, scgth = NA, scgnm = NA){
  bar <- fread(scbar, header = F)[[1]]
  fea <- fread(scfea, header = F)[[1]]
  library(stringr)
  fea <- str_to_upper(fea)
  mat <- readMM(scmat)
  mat <- as(mat, "CsparseMatrix")
  rownames(mat) <- fea
  colnames(mat) <- bar
  if(is.na(scgth)){
    sce <- SingleCellExperiment(assays = list(counts = mat),
                                rowData = DataFrame(gene_name = fea),
                                colData = DataFrame(barcodes = bar))
  }
  else{
    gth <- fread(scgth)
    gth <- gth[[scgnm]]
    sce <- SingleCellExperiment(assays = list(counts = mat),
                                rowData = DataFrame(gene_name = fea),
                                colData = DataFrame(barcodes = bar,
                                                    free_annotation = gth))
    colLabels(sce) <- gth
  }
  return(sce)
}

# Load scRNA-seq data by the users' definition: a txt file
Load_txtsc_to_SCE <- function(scmat, scgth = NA, scgnm = NA){
  mat <- fread(scmat)
  fea <- toupper(mat[[1]])
  mat <- as.matrix(mat[, -1])
  rownames(mat) <- fea
  bar <- c(paste("cell", 1:ncol(mat), sep = "_"))
  colnames(mat) <- bar
  mat <- as(mat, "CsparseMatrix")
  if(is.na(scgth)){
    sce <- SingleCellExperiment(assays = list(counts = mat),
                                rowData = DataFrame(gene_name = fea),
                                colData = DataFrame(barcodes = bar))
  }
  else{
    gth <- fread(scgth)
    gth <- gth[[scgnm]]
    sce <- SingleCellExperiment(assays = list(counts = mat),
                                rowData = DataFrame(gene_name = fea),
                                colData = DataFrame(barcodes = bar,
                                                    free_annotation = gth))
    colLabels(sce) <- gth
  }
  return(sce)
}


Load_h5adsc_to_SCE <- function(scmat, scgnm = NA){
  scmat <- h5read(scmat, '/')
  if("raw" %in% attr(scmat, "names")){
    X <- scmat$raw$X
    var <- scmat$raw$var
  } else{
    X <- scmat$X
    var <- scmat$var
  }
  # generate factors using categories
  for(name in attr(var[["__categories"]], "names")){
    if(length(var[[name]]) >= length(var[["__categories"]][[name]])){
      var[[name]] <- factor(var[[name]], labels = var[["__categories"]][[name]])
    }
  }
  var[["__categories"]] <- NULL
  obs <- scmat$obs
  for(name in attr(obs[["__categories"]], "names")){
    if(length(obs[[name]]) >= length(obs[["__categories"]][[name]]))
      obs[[name]] <- factor(obs[[name]], labels = obs[["__categories"]][[name]])
  }
  obs[["__categories"]] <- NULL
  dims <- c(length(scmat$var[["_index"]]), length(scmat$obs[["_index"]]))
  dat <- sparseMatrix(i = X$indices[] + 1,
                      p = X$indptr[],
                      x = as.numeric(X$data[]),
                      dims = dims,
                      dimnames = list(as.character(scmat$var[["_index"]]), as.character(scmat$obs[["_index"]])),
                      repr = "C")
  sce <- SingleCellExperiment(assays = list(counts = dat),
                              rowData = DataFrame(data.frame(var)),
                              colData = DataFrame(data.frame(obs)))

  # add ground truth
  if(!is.na(scgnm)){
    sce@colData$free_annotation <- obs[[scgnm]]
  }
  return(sce)
}
