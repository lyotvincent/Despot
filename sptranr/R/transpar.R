library(rhdf5)
library(Matrix)
library(png)
library(rjson)
library(data.table)

# Initialize a h5spt file at the beginning of pipline
Init_spt <- function(sptFile){
  # create a temp file
  h5createFile(sptFile)

  # create the metadata
  rhdf5::h5createGroup(sptFile, "metadata")

  # create first group `matrix`
  rhdf5::h5createGroup(sptFile, "matrix")
  # create `features`
  rhdf5::h5createGroup(sptFile, "matrix/features")
  # create `is_HVG`
  rhdf5::h5createGroup(sptFile, "matrix/features/is_HVG")
  # create `idents`
  rhdf5::h5createGroup(sptFile, "matrix/idents")
  # create `deconv`
  rhdf5::h5createGroup(sptFile, "matrix/deconv")
  # create `benchmark`
  rhdf5::h5createGroup(sptFile, "matrix/benchmark")

  # create second group `sptimages`
  rhdf5::h5createGroup(sptFile, "sptimages")
  # create `scalefactors`
  rhdf5::h5createGroup(sptFile, "sptimages/scalefactors")
  # create `coordinates`
  rhdf5::h5createGroup(sptFile, "sptimages/coordinates")

}

# Recording configs in sptFile
Save_cfg_to_spt <- function(sptFile, params){
  # create group `configs`
  rhdf5::h5createGroup(sptFile, "configs")
  for(par in names(params)){
    print(class(params[[par]]))
    Create_spt_array1d(sptFile, params[[par]], paste0("configs/", par), class(params[[par]]))
  }
  message("Configs are Loaded.")
}


# Save ground.truth to spt
Save_gt_to_spt <- function(sptFile, ground.truth, ground.sep=',', ground.name='ground_truth'){
  # Save ground_truth
  print(ground.truth)
  if(file.exists(ground.truth)){
    gtable <- fread(ground.truth, header = T)
    rownames(gtable) <- gtable[[1]]
    gtable <- gtable[order(rownames(gtable)),]
    gt <- gtable[[ground.name]]
    Create_spt_array1d(sptFile, gt, "matrix/idents/ground_truth", "character")
  }else{
    message("No ground truth for this data.")
  }
}

# Save 10X data to spt
Save_10X_to_spt <- function(dir, sptFile, filtered.matrix = T, name = "",
                            ground.truth = NULL, ground.sep = ",",
                            ground.name = "ground_truth", save.hires = F){
  library(stringr)
  # Save name
  if(length(name) != 0){
    Create_spt_array1d(sptFile, name, "name", "character")
  }

  # Save ground_truth
  if(!is.null(ground.truth)){
    Save_gt_to_spt(sptFile, ground.truth = ground.truth, ground.name = ground.name)
  }
  # whether .h5 file in the fold, if not load data from subdictionary: feature_bc_matrix
  h5dir <- ifelse(filtered.matrix,
                  "/filtered_feature_bc_matrix.h5",
                  "/raw_feature_bc_matrix.h5")
  h5dir <- paste0(dir, h5dir)
  if(file.exists(h5dir)){
    h5_obj <- rhdf5::h5read(h5dir, 'matrix')
  }
  else{
    h5dir <- ifelse(filtered.matrix,
                    "/filtered_feature_bc_matrix",
                    "/raw_feature_bc_matrix")
    h5dir <- paste0(dir, h5dir)
    if(file.exists(h5dir)){
      h5_obj <- list()
      h5_obj$barcodes <- fread(paste0(h5dir, "/barcodes.tsv.gz"), header = F)[[1]]
      h5_obj$features <- fread(paste0(h5dir, "/features.tsv.gz"), header = F)
      colnames(h5_obj$features) <- "name"
      mtx <- Matrix::readMM(paste0(h5dir, "/matrix.mtx.gz"))
      mtx <- as(mtx, "CsparseMatrix")
      mtx@Dimnames[[1]] <- toupper(h5_obj$features$names)
      mtx@Dimnames[[2]] <- toupper(h5_obj$barcodes)
      mtx <- mtx[, order(mtx@Dimnames[[2]])]
      h5_obj$barcodes <- toupper(h5_obj$barcodes[order(h5_obj$barcodes)])
      h5_obj$data <- mtx@x
      h5_obj$indices <- mtx@i
      h5_obj$indptr <- mtx@p
      h5_obj$shape <- mtx@Dim
      rm(mtx)
    }
    else{
      message("No supported files in the 10X folder, exit.")
      return(0)
    }
  }


  # save images
  library(Seurat)
  imgdir <- paste0(dir, "/spatial")
  # read in lowres, scalefactors, positions
  sf <- fromJSON(file = paste0(imgdir, "/scalefactors_json.json"))
  if(save.hires){
    loc <- H5Fopen(sptFile)
    if(H5Lexists(loc, "sptimages/hires")){
      H5Ldelete(loc, "sptimages/hires")
    }
    H5close()
    hires <- readPNG(paste0(imgdir, "/tissue_hires_image.png"))
    rhdf5::h5createDataset(sptFile, "sptimages/hires",
                           dim(hires), storage.mode = "double",
                           chunk = c(50, 50, 3))
    rhdf5::h5write(hires, sptFile, "sptimages/hires")
  }

  img <- Read10X_Image(imgdir, filter.matrix = F)
  # incorrect in Read10X_Image
  img@scale.factors$spot <- sf$spot_diameter_fullres
  image <- img@image

  loc <- H5Fopen(sptFile)
  if(H5Lexists(loc, "sptimages/lowres")){
    H5Ldelete(loc, "sptimages/lowres")
  }
  H5close()
  rhdf5::h5createDataset(sptFile, "sptimages/lowres",
                         dim(image), storage.mode = "double",
                         chunk = c(50, 50, 3))
  rhdf5::h5write(image, sptFile, "sptimages/lowres")

  coord <- img@coordinates
  coord <- coord[order(rownames(coord)), ]
  if(filtered.matrix){
    coord <- coord[rownames(coord) %in% h5_obj$barcodes, ]
    coord$tissue <- rep(1, length(h5_obj$barcodes))
  }

  # save names
  for(fac in attr(img@scale.factors, "names")){
    Create_spt_array1d(sptFile, img@scale.factors[[fac]],
                       paste0("sptimages/scalefactors/", fac),"double")
  }

  # save gene expression
  Create_spt_array1d(sptFile, toupper(h5_obj$barcodes), "matrix/barcodes", "character")
  Create_spt_array1d(sptFile, h5_obj$data, "matrix/data", "integer")
  Create_spt_array1d(sptFile, h5_obj$indices, "matrix/indices", "integer")
  Create_spt_array1d(sptFile, h5_obj$indptr, "matrix/indptr", "integer")
  Create_spt_array1d(sptFile, h5_obj$shape, "matrix/shape", "integer")

  for(fea in attr(h5_obj$features, "names")){
    Create_spt_array1d(sptFile, toupper(h5_obj$features[[fea]]),
                       paste0("matrix/features/", fea), "character")
  }

  # save spatial_locs
  for(coo in colnames(coord)){
    Create_spt_array1d(sptFile, coord[[coo]],
                       paste0("sptimages/coordinates/", coo), "integer")
  }
  index <- rownames(coord)
  Create_spt_array1d(sptFile, toupper(index),
                     "sptimages/coordinates/index", "character")
  Create_spt_array1d(sptFile, img@spot.radius, "sptimages/spot_radius", "double")
}

# Save SlideSeq data to spt
Save_SlideSeq_to_spt <- function(dir, sptFile){
  countdir <- paste0(dir, "/Puck.count")
  idxdir <- paste0(dir, "/Puck.idx")
  count <- read.table(countdir, header = T, sep = ',')
  idx <- read.table(idxdir, header = T, sep = ",")
  features <- data.frame(gene_id = count[,1], gene_name = count[,2])
  barcodes <- data.frame(barcodes = idx[,2])
  coord <- data.frame(x = idx[,3], y = idx[,4], row.names = idx[,2])
  count <- count[,-1:-2]
  rownames(count) <- features$gene_name
  colnames(count) <- barcodes$barcodes
  count <- Matrix(data = as.matrix(count))

  Create_spt_array1d(sptFile, barcodes$barcodes, "matrix/barcodes", "character")
  Create_spt_array1d(sptFile, count@x, "matrix/data", "integer")
  Create_spt_array1d(sptFile, count@i, "matrix/indices", "integer")
  Create_spt_array1d(sptFile, count@p, "matrix/indptr", "integer")
  Create_spt_array1d(sptFile, count@Dim, "matrix/shape", "integer")

  for(fea in attr(features, "names")){
    Create_spt_array1d(sptFile, features[[fea]],
                       paste0("matrix/features/", fea), "character")
  }

  for(coo in attr(coord, "names")){
    Create_spt_array1d(sptFile, coord[[coo]],
                       paste0("sptimages/coordinates/", coo), "integer")
  }
  index <- rownames(coord)
  Create_spt_array1d(sptFile, index,
                     "sptimages/coordinates/index", "character")
}



cal_chunk <- function(dimlen){
  chunk <- ifelse(dimlen < 320000, as.integer((dimlen + 3)/4), 80000)
  return(chunk)
}

Create_spt_group <- function(sptFile, group){
  loc <- H5Fopen(sptFile)
  if(!H5Lexists(loc, group)){
    h5createGroup(sptFile, sptloc)
  }
  H5close()
}

Create_spt_array1d <- function(sptFile, arr, sptloc, mode){
  arr_len <- length(arr)
  arr_chunk <- cal_chunk(arr_len)
  if(mode == "character"){
    sz <- max(na.omit(nchar(arr)))
    if(sz == 0){  # empty character
      sz <- 1
    }
  }
  loc <- H5Fopen(sptFile)
  if(H5Lexists(loc, sptloc)){
    H5Ldelete(loc, sptloc)
  }
  H5close()
  rhdf5::h5createDataset(sptFile, sptloc,
                         dims = arr_len,
                         storage.mode = mode,
                         size = sz,
                         chunk = arr_chunk)
  rhdf5::h5write(arr, sptFile, sptloc)
}

# Load spt to Seurat
Load_spt_to_Seurat <- function(sptFile, imgdir, h5data='matrix'){
  library(Seurat)
  seu_obj <- Load_h5_to_Seurat(sptFile, h5data)
  # h5img <- rhdf5::h5read(sptFile, "sptimages")
  seu_obj <- Load_img_to_Seurat(imgdir, seu_obj)
  return(seu_obj)
}

Load_h5_to_Seurat <- function(sptFile, h5data = 'matrix', assay = "Spatial"){
  library(Seurat)
  library(stringr)
  # find if need to do subset
  datas <- strsplit(h5data, '/')[[1]]
  data0 <- datas[1]
  h5mat <- rhdf5::h5read(sptFile, data0)
  feature.name <- as.character(toupper(h5mat$features$name))
  barcodes <- as.character(toupper(h5mat$barcodes))
  dat <- sparseMatrix(i = h5mat$indices[] + 1,
                      p = h5mat$indptr[],
                      x = as.numeric(h5mat$data[]),
                      dims = h5mat$shape[],
                      dimnames = list(feature.name, barcodes),
                      repr = "C")
  seu <- CreateSeuratObject(counts = dat, project = "Seurat", assay = assay)
  # for(fea in c('feature_type', 'genome', 'id', 'name')){
  #   seu[[assay]]@meta.features[fea] <- h5mat$features[fea]
  # }
  # if existing subset
  if(length(datas) > 1){
    subset <- datas[length(datas)]
    subarcodes <- h5mat[[subset]]$barcodes
    seu <- seu[,subarcodes]  # do subset
  }
  return(seu)
}

Load_img_to_Seurat <- function(imgdir, seu,
                               assay = "Spatial", imgname = "slice1"){
  library(Seurat)
  img <- Read10X_Image(imgdir)
  Seurat::DefaultAssay(object = img) <- assay
  img <- img[colnames(x = seu)]
  seu[[imgname]] <- img
  return(seu)
}

# Load data from 10X and convert to Seurat Object
Load10Xh5_to_Seurat <- function(dir, filtered.matrix = T,
                                assay = "Spatial", has.img = T,
                                imgdir = "spatial", imgname = "slice1"){
  library(Seurat)
  h5dir <- ifelse(filtered.matrix,
                  "/filtered_feature_bc_matrix.h5",
                  "/raw_feature_bc_matrix.h5")
  h5dir <- paste0(dir, h5dir)
  # call `Load_h5_to_Seurat` to read h5
  seu <- Load_h5_to_Seurat(h5dir, assay)

  # call `Load_img_to_Seurat` to read img
  if(has.img){
    imgdir <- paste0(dir, "/", imgdir)
    seu <- Load_img_to_Seurat(imgdir, seu, assay, imgname)
  }
  return(seu)
}


Save_spt_from_Seurat <- function(sptFile, seu, h5data = 'matrix'){

  sptloc = paste0(h5data, '/idents/Seurat')
  # create `idents` in matrix
  Create_spt_array1d(sptFile,
                     arr = as.numeric(seu$seurat_clusters),
                     sptloc = sptloc,
                     mode = "integer")
}


Save_spt_from_BayesSpace <-function (sptFile, Bayes,
                                     save.enhanced = F, h5data = 'matrix'){

  sptloc = paste0(h5data, '/idents/BayesSpace')
  # create `idents` in matrix
  Create_spt_array1d(sptFile,
                     arr = Bayes@colData@listData$spatial.cluster,
                     sptloc = sptloc,
                     mode = "integer")

  # create `is_HVG` group in features, which is differed by methods
  sptloc = paste0(h5data, '/features/is_HVG/BayesSpace')
  Create_spt_array1d(sptFile,
                     arr = Bayes@rowRanges@elementMetadata@listData$is.HVG,
                     sptloc = sptloc,
                     mode = "logical")

  # save enhance
  if(save.enhanced == TRUE){
    rhdf5::h5createGroup(sptFile, "sptimages/enhance")
    for(enh in attr(Bayes@colData@listData, "names")){
      Create_spt_array1d(sptFile,
                         arr = Bayes@colData@listData[[enh]],
                         sptloc = paste0("sptimages/enhance/", gsub("\\.", "_", enh)),
                         mode = "double")
    }
    Create_spt_array1d(sptFile,
                       arr = Bayes@colData@rownames,
                       sptloc = "sptimages/enhance/spotnames",
                       mode = "character")
  }
}

Save_spt_from_SPARK <- function(sptFile, sparkX, h5data = 'matrix'){

  rhdf5::h5createGroup(sptFile, paste0(h5data, "/features/is_HVG/SPARK"))

  Create_spt_array1d(sptFile,
                     arr = attr(sparkX$stats, "dimnames")[[1]],
                     sptloc = paste0(h5data, "/features/is_HVG/SPARK/name"),
                     mode = "character")

  Create_spt_array1d(sptFile,
                     arr = sparkX$res_mtest$combinedPval,
                     sptloc = paste0(h5data, "/features/is_HVG/SPARK/pval"),
                     mode = "double")
  Create_spt_array1d(sptFile,
                     arr = sparkX$res_mtest$adjustedPval,
                     sptloc = paste0(h5data, "/features/is_HVG/SPARK/qval"),
                     mode = "double")
}


Save_spt_from_Seurat.markers <- function(sptFile, de_markers, h5data = 'matrix'){
  rhdf5::h5createGroup(sptFile, paste0(h5data, "/features/is_HVG/Seurat"))
  Create_spt_array1d(sptFile,
                     arr = de_markers$gene,
                     sptloc = paste0(h5data, "/features/is_HVG/Seurat/name"),
                     mode = "character")
  Create_spt_array1d(sptFile,
                     arr = as.numeric(de_markers$cluster),
                     sptloc = paste0(h5data, "/features/is_HVG/Seurat/cluster"),
                     mode = "integer")
  Create_spt_array1d(sptFile,
                     arr = de_markers$p_val,
                     sptloc = paste0(h5data, "/features/is_HVG/Seurat/pval"),
                     mode = "double")
  Create_spt_array1d(sptFile,
                     arr = de_markers$p_val_adj,
                     sptloc = paste0(h5data, "/features/is_HVG/Seurat/qval"),
                     mode = "double")
}


Save_spt_from_Giotto_clu <- function(sptFile, gobj, h5data = 'matrix'){
  sptloc = paste0(h5data, '/idents/Giotto')
  # create `idents` in matrix
  Create_spt_array1d(sptFile,
                     arr = as.numeric(gobj@cell_metadata$leiden_clus),
                     sptloc = sptloc,
                     mode = "integer")
}

Save_spt_from_Giotto_est <- function(sptFile, spatialgenes, h5data = 'matrix'){

  rhdf5::h5createGroup(sptFile, paste0(h5data, "/features/is_HVG/Giotto"))

  Create_spt_array1d(sptFile,
                     arr = spatialgenes$genes,
                     sptloc = paste0(h5data, "/features/is_HVG/Giotto/name"),
                     mode = "character")

  Create_spt_array1d(sptFile,
                     arr = spatialgenes$p.value,
                     sptloc = paste0(h5data, "/features/is_HVG/Giotto/pval"),
                     mode = "double")
  Create_spt_array1d(sptFile,
                     arr = spatialgenes$adj.p.value,
                     sptloc = paste0(h5data, "/features/is_HVG/Giotto/qval"),
                     mode = "double")
}

Save_spt_from_Giotto_dcv <- function(sptFile, gobj, h5data = 'matrix'){
  result <- gobj@spatial_enrichment$DWLS

  barcodes <- result[[1]]
  result <- result[, -1]
  cell_type <- colnames(result)
  result <- as.matrix(result)
  result <- as(result, "dgeMatrix")
  h5createGroup(sptFile, paste0(h5data, '/deconv/Giotto'))

  # save 1d weights
  Create_spt_array1d(sptFile,
                     arr = result@x,
                     sptloc = paste0(h5data, '/deconv/Giotto/weights'),
                     mode = 'double')
  # save shape
  Create_spt_array1d(sptFile,
                     arr = result@Dim,
                     sptloc = paste0(h5data, '/deconv/Giotto/shape'),
                     mode = 'integer')
  # save dim names
  Create_spt_array1d(sptFile,
                     arr = barcodes,
                     sptloc = paste0(h5data, '/deconv/Giotto/barcodes'),
                     mode = 'character')
  Create_spt_array1d(sptFile,
                     arr = cell_type,
                     sptloc = paste0(h5data, '/deconv/Giotto/cell_type'),
                     mode = 'character')

}

# spt to SingleCellExperiment object
Load_spt_to_SCE <- function(sptFile, h5data = 'matrix'){
  library(SingleCellExperiment)
  library(SummarizedExperiment)

  # find if need to do subset
  datas <- strsplit(h5data, '/')[[1]]
  data0 <- datas[1]

  h5mat <- rhdf5::h5read(sptFile, data0)
  feature.name <- as.character(toupper(h5mat$features$name))
  barcodes <- as.character(toupper(h5mat$barcodes))
  h5img <- rhdf5::h5read(sptFile, "sptimages")
  dat <- sparseMatrix(i = h5mat$indices[] + 1,
                      p = h5mat$indptr[],
                      x = as.numeric(h5mat$data[]),
                      dims = h5mat$shape[],
                      dimnames = list(feature.name, barcodes),
                      repr = "C")
  colData <- DataFrame(spot = h5img$coordinates$index,
                     in_tissue = h5img$coordinates$tissue,
                     row = h5img$coordinates$row,
                     col = h5img$coordinates$col,
                     imagerow = h5img$coordinates$imagerow,
                     imagecol = h5img$coordinates$imagecol)
  colData <- colData[order(colData$spot),]
  colData <- colData[colData$spot %in% h5mat$barcodes, ]
  for(ident in attr(h5mat$idents, "names")){
    id <- data.frame(h5mat$idents[[ident]])
    colnames(id) <- ident
    colData <- cbind(colData, id)
  }
  # some datasets may not provide gene_id.
  if("id" %in% attr(h5mat$features, "names")){
    rowData <- DataFrame(gene_id = h5mat$features$id,
                        gene_name = h5mat$features$name)
  }
  else{
    rowData <- DataFrame(gene_name = h5mat$features$name)
  }
  sce <- SingleCellExperiment(assays = list(counts = dat),
                             rowData = rowData,
                             colData = colData)

  # if existing subset
  if(length(datas) > 1){
    subset <- datas[length(datas)]
    subarcodes <- h5mat[[subset]]$barcodes
    sce <- sce[,subarcodes]  # do subset
  }
  return(sce)
}


# spt to Giotto object
Load_spt_to_Giotto <- function(sptFile, h5data = 'matrix',
                               python_path = '~/anaconda3/envs/spatial/bin/python',
                               temp_dir = 'h5ads'){
  library(Giotto)
  library(SeuratObject)
  library(stringr)
  ginst <- createGiottoInstructions(python_path = python_path,
                                    save_dir = temp_dir,
                                    save_plot = FALSE,
                                    show_plot = FALSE)

  # find if need to do subset
  datas <- strsplit(h5data, '/')[[1]]
  data0 <- datas[1]

  # read the origin file with features, barcodes and images
  h5_obj <- rhdf5::h5read(sptFile, data0)
  feature.name <- toupper(h5_obj$features$name)
  barcodes <- toupper(h5_obj$barcodes)
  h5img <- rhdf5::h5read(sptFile, "sptimages")

  # load spatial locations
  spatial_locs <- data.frame(row.names = h5img$coordinates$index,
                        sdimx = h5img$coordinates$imagerow,
                        sdimy = h5img$coordinates$imagecol)
  spatial_locs <- spatial_locs[order(rownames(spatial_locs)),]

  # load expression matrix as dgCMatrix
  dat <- sparseMatrix(i = h5_obj$indices[] + 1,
                      p = h5_obj$indptr[],
                      x = as.numeric(h5_obj$data[]),
                      dims = h5_obj$shape[],
                      dimnames = list(feature.name, barcodes),
                      repr = "C")
  obj <- CreateSeuratObject(dat)

  # if existing subset
  if(length(datas) > 1){
    subset <- datas[length(datas)]
    subarcodes <- h5_obj[[subset]]$barcodes
    obj <- obj[,subarcodes]  # do subset using Seurat
    spatial_locs <- spatial_locs[subarcodes,]  # do subset for spatial_locs
  }
  gobj <- createGiottoObject(raw_exprs = obj@assays$RNA@counts,
                             spatial_locs = spatial_locs,
                             instructions = ginst)
  return(gobj)
}

# spt to SpatialExperiment object
Load_spt_to_SE <- function(sptFile, h5data = 'matrix'){

  library(SpatialExperiment)
  # find if need to do subset
  datas <- strsplit(h5data, '/')[[1]]
  data0 <- datas[1]

  h5mat <- rhdf5::h5read(sptFile, data0)
  feature.name <- as.character(toupper(h5mat$features$name))
  barcodes <- as.character(toupper(h5mat$barcodes))
  h5img <- rhdf5::h5read(sptFile, "sptimages")
  # read in spatial coordinates
  coord <- data.frame(h5img$coordinates,
                      row.names = h5img$coordinates$index)
  coord <- coord[order(h5img$coordinates$index),]
  # read in counts
  dat <- sparseMatrix(i = h5mat$indices[] + 1,
                      p = h5mat$indptr[],
                      x = as.numeric(h5mat$data[]),
                      dims = h5mat$shape[],
                      dimnames = list(feature.name, barcodes),
                      repr = "C")
  dat_assay <- S4Vectors::SimpleList(dat)
  # construct observation & feature metadata
  rowData <- DataFrame(symbol = feature.name)

  # barcodes data
  colData <- DataFrame(Barcode = barcodes)
  for(ident in attr(h5mat$idents, "names")){
    id <- data.frame(h5mat$idents[[ident]])
    colnames(id) <- ident
    colData <- cbind(colData, id)
  }
  # read in images
  rgb <- as.raster(h5img$lowres)
  rgb <- SpatialExperiment::SpatialImage(rgb)
  img <- S4Vectors::DataFrame(sample_id = "foo",
                              image_id = "lowres")
  img[["data"]] <- list(rgb)
  img[["scaleFactor"]] <- h5img$scalefactors$lowres
  spe <- SpatialExperiment(assays = list(counts = dat),
                           colData = colData, rowData = rowData, imgData = img,
                           spatialData = DataFrame(coord),
                           sample_id = "foo")
  # if existing subset
  if(length(datas) > 1){
    subset <- datas[length(datas)]
    subarcodes <- h5_obj[[subset]]$barcodes
    spe <- spe[,subarcodes]  # do subset
  }
  return(spe)
}

# spt to SPARK object
Load_spt_to_SPARK <- function(sptFile, h5data = 'matrix'){
  library(SPARK)

  # find if need to do subset
  datas <- strsplit(h5data, '/')[[1]]
  data0 <- datas[1]

  h5mat <- rhdf5::h5read(sptFile, data0)
  feature.name <- as.character(toupper(h5mat$features$name))
  barcodes <- as.character(toupper(h5mat$barcodes))

  dat <- sparseMatrix(i = h5mat$indices[] + 1,
                      p = h5mat$indptr[],
                      x = as.numeric(h5mat$data[]),
                      dims = h5mat$shape[],
                      dimnames = list(feature.name, barcodes),
                      repr = "C")
  h5img <- rhdf5::h5read(sptFile, "sptimages")
  spatial_locs <- data.frame(x = h5img$coordinates$row, y = h5img$coordinates$col)
  rownames(spatial_locs) <- colnames(dat)

  obj <- CreateSeuratObject(dat)

  # if existing subset
  if(length(datas) > 1){
    subset <- datas[length(datas)]
    subarcodes <- h5_obj[[subset]]$barcodes
    obj <- obj[,subarcodes]  # do subset using Seurat
    spatial_locs <- spatial_locs[subarcodes,]  # do subset for spatial_locs
  }
  spark <- CreateSPARKObject(counts = obj@assays$RNA@counts,
                             location = spatial_locs,
                             project = "SPARK",
                             percentage = 0.1,
                             min_total_counts = 10)

  ## total counts for each cell/spot
  spark@lib_size <- apply(spark@counts, 2, sum)
  return(spark)
}

# spt to SpatialRNA object
Load_spt_to_SpatialRNA <- function(sptFile){
  library(spacexr)
  library(SeuratObject)
  # find if need to do subset
  datas <- strsplit(h5data, '/')[[1]]
  data0 <- datas[1]

  h5mat <- rhdf5::h5read(sptFile, data0)
  feature.name <- as.character(toupper(h5mat$features$name))
  barcodes <- as.character(toupper(h5mat$barcodes))
  h5img <- rhdf5::h5read(sptFile, "sptimages")
  # create coords for SpatialRNA
  spatial_locs <- data.frame(x = h5img$coordinates$row,
                      y = h5img$coordinates$col,
                      row.names = barcodes)

  # genereate CsparseMatrix in counts for SpatialRNA
  dat <- sparseMatrix(i = h5mat$indices[] + 1,
                       p = h5mat$indptr[],
                       x = as.numeric(h5mat$data[]),
                       dims = h5mat$shape[],
                       dimnames = list(feature.name, barcodes),
                       repr = "C")

  obj <- CreateSeuratObject(dat, min.cells = 10)

  # if existing subset
  if(length(datas) > 1){
    subset <- datas[length(datas)]
    subarcodes <- h5_obj[[subset]]$barcodes
    obj <- obj[,subarcodes]  # do subset using Seurat
    spatial_locs <- spatial_locs[subarcodes,]  # do subset for spatial_locs
  }

  sr <- SpatialRNA(coords = spatial_locs,
                   counts = obj@assays$RNA@counts,
                   use_fake_coords = FALSE,
                   require_int = TRUE)
  return(sr)
}

# some sofrwares only support .tsv files, so we provide the interface
Save_spt_to_tsv <- function(sptFile, outdir){
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  h5mat <- rhdf5::h5read(sptFile, '/matrix')

  # genereate CsparseMatrix in counts
  counts <- sparseMatrix(i = h5mat$indices[] + 1,
                         p = h5mat$indptr[],
                         x = as.numeric(h5mat$data[]),
                         dims = h5mat$shape[],
                         dimnames = list(h5mat$features$name, h5mat$barcodes),
                         repr = "C")
  spmat <- as.matrix(counts)
  write.table(spmat, file = paste0(outdir, '/spcnt.tsv'), sep = '\t', quote = F)
}
