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

  # create second group `images`
  rhdf5::h5createGroup(sptFile, "sptimages")
  # create `scalefactors`
  rhdf5::h5createGroup(sptFile, "sptimages/scalefactors")
  # create `coordinates`
  rhdf5::h5createGroup(sptFile, "sptimages/coordinates")

}

# Save ground.truth to spt
Save_gt_to_spt <- function(sptFile, ground.truth, ground.sep=',', ground.name='ground_truth'){
  # Save ground_truth
  if(dir.exists(ground.truth)){
    message("No ground truth for this data.")
  }
  else{
    if(file.exists(ground.truth)){
      gtable <- fread(ground.truth, header = T)
      gt <- gtable[[ground.name]]
      Create_spt_array1d(sptFile, gt, "matrix/idents/ground_truth", "character")
    }else{
      message("No ground truth for this data.")
    }
  }
}

# Save 10X data to spt
Save_10X_to_spt <- function(dir, sptFile, filtered.matrix = T, name = "",
                            ground.truth = NULL, ground.sep = ",",
                            ground.name = "ground_truth", save.hires = F){
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
                    "/filtered_feature_bc_matrix/",
                    "/raw_feature_bc_matrix/")
    h5dir <- paste0(dir, h5dir)
    if(file.exists(h5dir)){
      h5_obj <- list()
      h5_obj$barcodes <- fread(paste0(h5dir, "barcodes.tsv.gz"), header = F)[[1]]
      h5_obj$features <- fread(paste0(h5dir, "features.tsv.gz"), header = F)
      colnames(h5_obj$features) <- "names"
      mtx <- Matrix::readMM(paste0(h5dir, "matrix.mtx.gz"))
      mtx <- as(mtx, "dgCMatrix")
      mtx@Dimnames[[1]] <- h5_obj$features$names
      mtx@Dimnames[[2]] <- h5_obj$barcodes
      mtx <- mtx[, order(mtx@Dimnames[[2]])]
      h5_obj$barcodes <- h5_obj$barcodes[order(h5_obj$barcodes)]
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
  Create_spt_array1d(sptFile, h5_obj$barcodes, "matrix/barcodes", "character")
  Create_spt_array1d(sptFile, h5_obj$data, "matrix/data", "integer")
  Create_spt_array1d(sptFile, h5_obj$indices, "matrix/indices", "integer")
  Create_spt_array1d(sptFile, h5_obj$indptr, "matrix/indptr", "integer")
  Create_spt_array1d(sptFile, h5_obj$shape, "matrix/shape", "integer")

  for(fea in attr(h5_obj$features, "names")){
    Create_spt_array1d(sptFile, h5_obj$features[[fea]],
                       paste0("matrix/features/", fea), "character")
  }

  # save spatial_locs
  for(coo in colnames(coord)){
    Create_spt_array1d(sptFile, coord[[coo]],
                       paste0("sptimages/coordinates/", coo), "integer")
  }
  index <- rownames(coord)
  Create_spt_array1d(sptFile, index,
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

Load_h5_to_Seurat <- function(h5dir, h5data = 'matrix', assay = "Spatial"){
  library(Seurat)
  if(h5data == 'matrix'){
    h5_obj <- rhdf5::h5read(h5dir, h5data)
    dat <- sparseMatrix(i = h5_obj$indices[] + 1,
                        p = h5_obj$indptr[],
                        x = as.numeric(h5_obj$data[]),
                        dims = h5_obj$shape[],
                        dimnames = list(h5_obj$features$id, h5_obj$barcodes),
                        repr = "C")
    seu <- CreateSeuratObject(counts = dat, project = "Seurat", assay = assay)
    for(fea in c('feature_type', 'genome', 'id', 'name')){
      seu[[assay]]@meta.features[fea] <- h5_obj$features[fea]
    }
  }
  else if(h5data == 'SpotClean_mat'){
    h5_obj <- rhdf5::h5read(h5dir, h5data)
    dat <- sparseMatrix(i = h5_obj$indices[] + 1,
                        p = h5_obj$indptr[],
                        x = as.numeric(h5_obj$data[]),
                        dims = h5_obj$shape[],
                        dimnames = list(h5_obj$features$name, h5_obj$barcodes),
                        repr = "C")
    seu <- CreateSeuratObject(counts = dat, project = "Seurat", assay = assay)
  }
  else if(h5data == 'SPCS_mat'){
    h5_obj <- rhdf5::h5read(h5dir, h5data)
    dat <- sparseMatrix(i = h5_obj$indices[] + 1,
                        p = h5_obj$indptr[],
                        x = as.numeric(h5_obj$data[]),
                        dims = h5_obj$shape[],
                        dimnames = list(h5_obj$features$name, h5_obj$barcodes),
                        repr = "C")
    seu <- CreateSeuratObject(counts = dat, project = "Seurat", assay = assay)
  }
  else if(h5data == 'SPCS_mat'){
    h5_obj <- rhdf5::h5read(h5dir, h5data)
    dat <- sparseMatrix(i = h5_obj$indices[] + 1,
                        p = h5_obj$indptr[],
                        x = as.numeric(h5_obj$data[]),
                        dims = h5_obj$shape[],
                        dimnames = list(h5_obj$feature_name, h5_obj$barcodes),
                        repr = "C")
    seu <- CreateSeuratObject(counts = dat, project = "Seurat", assay = assay)
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

  # library(Seurat)
  # assay <- "Spatial"
  #
  # # create a temp file
  # h5createFile(h5File)
  #
  # # create first group `matrix`
  # rhdf5::h5createGroup(h5File, "matrix")
  #
  # # create `barcodes` in `matrix`, from active.ident in `seurat_obj`
  # Create_spt_array1d(h5File,
  #                    arr = attr(seu_obj@active.ident, "names"),
  #                    sptloc = "matrix/barcodes",
  #                    mode = "character")
  #
  # # create `data` in `matrix`
  # data <- seu_obj@assays$Spatial@counts@x
  # dat_len <- length(data)
  # dat_chunk <- cal_chunk(dat_len)
  # rhdf5::h5createDataset(h5File, "matrix/data",
  #                        dims = dat_len,
  #                        storage.mode = "integer",
  #                        chunk = dat_chunk)
  # rhdf5::h5write(data, h5File, "matrix/data")
  #
  # # create second group `features` and `_all_tag_keys`
  # rhdf5::h5createGroup(h5File, "matrix/features")
  # rhdf5::h5createDataset(h5File, "matrix/features/_all_tag_keys",
  #                        dims = 1, storage.mode = "character", size = 6)
  # rhdf5::h5write("genome", h5File, "matrix/features/_all_tag_keys")
  #
  # # create `features`
  # feature <- seu_obj[[assay]]@meta.features
  # fea_len <- seu_obj[[assay]]@data@Dim[1]
  # fea_chunk <- cal_chunk(fea_len)
  # fea_loc <- "matrix/features/"
  # for(fea in colnames(feature)){
  #   flist <- c(t(feature[fea]))
  #   rhdf5::h5createDataset(h5File, paste0(fea_loc, fea),
  #                          dims = fea_len,
  #                          storage.mode = "character",
  #                          size = max(nchar(flist)),
  #                          chunk = fea_chunk)
  #   rhdf5::h5write(flist, h5File, paste0(fea_loc, fea))
  # }
  #
  # # create indices
  # ind <- seu_obj[['Spatial']]@counts@i
  # ind_len <- length(ind)
  # ind_chunk <- cal_chunk(ind_len)
  # rhdf5::h5createDataset(h5File, "matrix/indices",
  #                        dims = ind_len,
  #                        storage.mode = "integer",
  #                        chunk = ind_chunk)
  # rhdf5::h5write(ind, h5File, "matrix/indices")
  #
  # # create indptr
  # ptr <- seu_obj[["Spatial"]]@counts@p
  # ptr_len <- length(ptr)
  # ptr_chunk <- cal_chunk(ptr_len)
  # rhdf5::h5createDataset(h5File, "matrix/indptr",
  #                        dims = ptr_len,
  #                        storage.mode = "integer",
  #                        chunk = ptr_chunk)
  # rhdf5::h5write(ptr, h5File, "matrix/indptr")
  #
  # # create shape
  # shp <- seu_obj[['Spatial']]@counts@Dim
  # rhdf5::h5createDataset(h5File, "matrix/shape",
  #                        2, storage.mode = "integer")
  # rhdf5::h5write(shp, h5File, "matrix/shape")

  # # create images
  # rhdf5::h5createGroup(h5File, "images")
  # # find all slices
  # slices <- attr(seu_obj@images, "names")
  # for(sl in slices){
  #   slice <- seu_obj@images[[sl]]
  #   sl_path <- paste0("images/", sl)
  #   # create a group for single slice
  #   rhdf5::h5createGroup(h5File, sl_path)
  #
  #   # save image under the slice
  #   img <- slice@image
  #   img_path <- paste0(sl_path, "/image")
  #   img_chunk <- cal_chunk(length(img))
  #   rhdf5::h5createDataset(h5File, img_path,
  #                          dim(img), storage.mode = "double",
  #                          chunk = c(50, 50, 3))
  #   rhdf5::h5write(img, h5File, img_path)
  #
  #   # save scale factors under the slice
  #   scl_path <- paste0(sl_path, "/scalefactors")
  #   rhdf5::h5createGroup(h5File, scl_path)
  #   for(fac in attr(slice@scale.factors, "names")){
  #     fac_path <- paste0(scl_path, "/", fac)
  #     rhdf5::h5createDataset(h5File, fac_path,
  #                            1, storage.mode = "double")
  #     rhdf5::h5write(slice@scale.factors[[fac]], h5File, fac_path)
  #   }
  #
  #   # save coordinates under the slice
  #   coo_path <- paste0(sl_path, "/coordinates")
  #   rhdf5::h5createGroup(h5File, coo_path)
  #   for(coo in attr(slice@coordinates, "names")){
  #     pos <- slice@coordinates[[coo]]
  #     pos_path <- paste0(coo_path, "/", coo)
  #     rhdf5::h5createDataset(h5File, pos_path,
  #                            length(pos), storage.mode = "integer")
  #     rhdf5::h5write(pos, h5File, pos_path)
  #   }
  #
  #   # save spot radius
  #   rad_path <- paste0(sl_path, "/spot_radius")
  #   rhdf5::h5createDataset(h5File, rad_path,
  #                          1, storage.mode = "double")
  #   rhdf5::h5write(slice@spot.radius, h5File, rad_path)
  # }
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

# Seurat to SingleCellExperiment object
Tran_Seurat_to_SCE <-function(seu_obj){
  Bayes <- as.SingleCellExperiment(seu_obj)
  Bayes@rowRanges@elementMetadata@listData[["gene_id"]] <- seu_obj[['Spatial']]@meta.features$id
  Bayes@rowRanges@elementMetadata@listData[["gene_name"]] <- seu_obj[['Spatial']]@meta.features$name
  Bayes@colData@listData[["spot"]] <- attr(seu_obj@active.ident, "names")
  Bayes@colData@listData[["in_tissue"]] <- seu_obj@images$slice1@coordinates$tissue
  Bayes@colData@listData[["row"]] <- seu_obj@images[["slice1"]]@coordinates[["row"]]
  Bayes@colData@listData[["col"]] <- seu_obj@images[["slice1"]]@coordinates[["col"]]
  Bayes@colData@listData[["imagerow"]] <- seu_obj@images[["slice1"]]@coordinates[["imagerow"]]
  Bayes@colData@listData[["imagecol"]] <- seu_obj@images[["slice1"]]@coordinates[["imagecol"]]
  return(Bayes)
}

# spt to SingleCellExperiment object
Load_spt_to_SCE <- function(sptFile, h5data = 'matrix'){
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  if(h5data == 'matrix'){
    h5_obj <- rhdf5::h5read(sptFile, h5data)
    h5img <- rhdf5::h5read(sptFile, "sptimages")
    dat <- sparseMatrix(i = h5_obj$indices[] + 1,
                        p = h5_obj$indptr[],
                        x = as.numeric(h5_obj$data[]),
                        dims = h5_obj$shape[],
                        dimnames = list(NULL, NULL),
                        repr = "C")
    colData <- DataFrame(spot = h5img$coordinates$index,
                       in_tissue = h5img$coordinates$tissue,
                       row = h5img$coordinates$row,
                       col = h5img$coordinates$col,
                       imagerow = h5img$coordinates$imagerow,
                       imagecol = h5img$coordinates$imagecol)
    colData <- colData[order(colData$spot),]
    colData <- colData[colData$spot %in% h5_obj$barcodes, ]
    # some datasets may not provide gene_id.
    if("id" %in% attr(h5_obj$features, "names")){
      rowData <- DataFrame(gene_id = h5_obj$features$id,
                          gene_name = h5_obj$features$name)
    }
    else{
      rowData <- DataFrame(gene_name = h5_obj$features$name)
    }
    se <- SummarizedExperiment(assays = list(counts = dat),
                               rowData = rowData,
                               colData = colData)
    sce <- as(se, "SingleCellExperiment")
    colnames(sce) <- h5_obj$barcodes
    rownames(sce) <- h5_obj$features$name
  }
  else if(h5data == "SpotClean_mat"){
    h5_obj <- rhdf5::h5read(sptFile, h5data)
    h5img <- rhdf5::h5read(sptFile, "sptimages")
    dat <- sparseMatrix(i = h5_obj$indices[] + 1,
                        p = h5_obj$indptr[],
                        x = as.numeric(h5_obj$data[]),
                        dims = h5_obj$shape[],
                        dimnames = list(NULL, NULL),
                        repr = "C")
    colData <- DataFrame(spot = h5img$coordinates$index,
                         in_tissue = h5img$coordinates$tissue,
                         row = h5img$coordinates$row,
                         col = h5img$coordinates$col,
                         imagerow = h5img$coordinates$imagerow,
                         imagecol = h5img$coordinates$imagecol)
    colData <- colData[order(colData$spot),]
    colData <- colData[colData$spot %in% h5_obj$barcodes, ]
    se <- SummarizedExperiment(assays = list(counts = dat),
                               rowData = DataFrame(gene_name = h5_obj$features$name),
                               colData = colData)
    sce <- as(se, "SingleCellExperiment")
    colnames(sce) <- h5_obj$barcodes
    rownames(sce) <- h5_obj$features$name
  }
  else if(h5data == "SPCS_mat"){
    h5_obj <- rhdf5::h5read(sptFile, h5data)
    h5img <- rhdf5::h5read(sptFile, "sptimages")
    dat <- sparseMatrix(i = h5_obj$indices[] + 1,
                        p = h5_obj$indptr[],
                        x = as.numeric(h5_obj$data[]),
                        dims = h5_obj$shape[],
                        dimnames = list(NULL, NULL),
                        repr = "C")
    colData <- DataFrame(spot = h5img$coordinates$index,
                         in_tissue = h5img$coordinates$tissue,
                         row = h5img$coordinates$row,
                         col = h5img$coordinates$col,
                         imagerow = h5img$coordinates$imagerow,
                         imagecol = h5img$coordinates$imagecol)
    colData <- colData[order(colData$spot),]
    colData <- colData[colData$spot %in% h5_obj$barcodes, ]
    se <- SummarizedExperiment(assays = list(counts = dat),
                               rowData = DataFrame(gene_name = h5_obj$features$name),
                               colData = colData)
    sce <- as(se, "SingleCellExperiment")
    colnames(sce) <- h5_obj$barcodes
    rownames(sce) <- h5_obj$features$name
  }
  return(sce)
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


# spt to Giotto object
Load_spt_to_Giotto <- function(sptFile, h5data = 'matrix',
                               python_path = '~/anaconda3/envs/spatial/bin/python',
                               temp_dir = 'h5ads'){
  library(Giotto)
  ginst <- createGiottoInstructions(python_path = python_path,
                                    save_dir = temp_dir,
                                    save_plot = FALSE,
                                    show_plot = FALSE)
  h5_obj <- rhdf5::h5read(sptFile, h5data)
  h5img <- rhdf5::h5read(sptFile, "sptimages")
  colData <- data.frame(row.names = h5img$coordinates$index,
                        sdimx = h5img$coordinates$imagerow,
                        sdimy = h5img$coordinates$imagecol)
  colData <- colData[order(rownames(colData)),]
  if(h5data == 'matrix'){
    dat <- sparseMatrix(i = h5_obj$indices[] + 1,
                        p = h5_obj$indptr[],
                        x = as.numeric(h5_obj$data[]),
                        dims = h5_obj$shape[],
                        dimnames = list(h5_obj$features$id, h5_obj$barcodes),
                        repr = "C")
    gobj <- createGiottoObject(raw_exprs = dat,
                               spatial_locs = colData,
                               instructions = ginst)
  }
  else if(h5data == "SpotClean_mat"){
    dat <- sparseMatrix(i = h5_obj$indices[] + 1,
                        p = h5_obj$indptr[],
                        x = as.numeric(h5_obj$data[]),
                        dims = h5_obj$shape[],
                        dimnames = list(h5_obj$features$name, h5_obj$barcodes),
                        repr = "C")
    gobj <- createGiottoObject(raw_exprs = dat,
                               spatial_locs = colData,
                               instructions = ginst)
  }
  else if(h5data == "SPCS_mat"){
    dat <- sparseMatrix(i = h5_obj$indices[] + 1,
                        p = h5_obj$indptr[],
                        x = as.numeric(h5_obj$data[]),
                        dims = h5_obj$shape[],
                        dimnames = list(h5_obj$features$name, h5_obj$barcodes),
                        repr = "C")
    gobj <- createGiottoObject(raw_exprs = dat,
                               spatial_locs = colData,
                               instructions = ginst)
  }
  return(gobj)
}

# spt to SpatialExperiment object
Load_spt_to_SE <- function(sptFile, h5data = 'matrix'){
  library(SpatialExperiment)
  h5mat <- rhdf5::h5read(sptFile, h5data)
  h5img <- rhdf5::h5read(sptFile, "sptimages")
  # read in spatial coordinates
  coord <- data.frame(h5img$coordinates,
                      row.names = h5mat$barcodes)
  # read in counts
  dat <- sparseMatrix(i = h5mat$indices[] + 1,
                      p = h5mat$indptr[],
                      x = as.numeric(h5mat$data[]),
                      dims = h5mat$shape[],
                      dimnames = list(NULL, NULL),
                      repr = "C")
  dat_assay <- S4Vectors::SimpleList(dat)
  # construct observation & feature metadata
  if(h5data == 'matrix'){
    rowData <- DataFrame(symbol = h5mat$features$name,
                    id = h5mat$features$id)
  }
  else{
    rowData <- DataFrame(symbol = h5mat$features$name)
  }

  # barcodes data
  colData <- DataFrame(Barcode = h5mat$barcodes)

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
  return(spe)
}

# spt to SPARK object
Load_spt_to_SPARK <- function(sptFile){
  library(SPARK)
  h5_obj <- rhdf5::h5read(sptFile, '/')
  h5mat <- h5_obj$matrix
  sdat <- sparseMatrix(i = h5mat$indices[] + 1,
                      p = h5mat$indptr[],
                      x = as.numeric(h5mat$data[]),
                      dims = h5mat$shape[],
                      dimnames = list(NULL, NULL),
                      repr = "C")
  dat <- as.data.frame(sdat)
  colnames(dat) <- h5mat$barcodes
  rownames(dat) <- h5mat$features$id
  h5img <- h5_obj$sptimages
  loc <- data.frame(x = h5img$coordinates$row, y = h5img$coordinates$col)
  rownames(loc) <- colnames(dat)
  spark <- CreateSPARKObject(counts = dat,
                             location = loc,
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
  h5_obj <- rhdf5::h5read(sptFile, '/')
  h5mat <- h5_obj$matrix

  # genereate dgCMatrix in counts for SpatialRNA
  counts <- sparseMatrix(i = h5mat$indices[] + 1,
                       p = h5mat$indptr[],
                       x = as.numeric(h5mat$data[]),
                       dims = h5mat$shape[],
                       dimnames = list(h5mat$features$name, h5mat$barcodes),
                       repr = "C")

  assay <- CreateAssayObject(counts, min.cells = 10)
  # create nUMI for SpatialRNA
  nUMI <- apply(assay@counts, 2, sum)

  h5img <- h5_obj$sptimages
  # create coords for SpatialRNA
  coord <- data.frame(x = h5img$coordinates$row,
                      y = h5img$coordinates$col)
  rownames(coord) <- h5mat$barcodes

  sr <- SpatialRNA(coords = coord,
                   counts = assay@counts,
                   nUMI = nUMI,
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

  # genereate dgCMatrix in counts
  counts <- sparseMatrix(i = h5mat$indices[] + 1,
                         p = h5mat$indptr[],
                         x = as.numeric(h5mat$data[]),
                         dims = h5mat$shape[],
                         dimnames = list(h5mat$features$name, h5mat$barcodes),
                         repr = "C")
  spmat <- as.matrix(counts)
  write.table(spmat, file = paste0(outdir, '/spcnt.tsv'), sep = '\t', quote = F)
}
