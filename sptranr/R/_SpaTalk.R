library(SpaTalk)
source("sptranr/R/transpar.R")

Deconv_SpaTalk <- function(spe, sce, species = "Mouse", is_cell = F){
  # create SpaTalk Object from SpatialExperiment and SingleCellExperiment object
  st_data <- spe@assays@data@listData$counts
  colnames(st_data) <- spe@colData$Barcode
  rownames(st_data) <- spe@rowRanges@elementMetadata$symbol
  sc_data <- sce@assays@data@listData$counts
  rownames(sc_data) <- rownames(sce)
  colnames(sc_data) <- colnames(sce)
  sc_celltype <- sce$free_annotation
  if(is_cell){
    st_meta <- data.frame(cell = spe@colData$Barcode,
                          x = spe@int_colData@listData$spatialData$row,
                          y = spe@int_colData@listData$spatialData$col)
  }else{
    st_meta <- data.frame(spot = spe@colData$Barcode,
                          x = spe@int_colData@listData$spatialData$row,
                          y = spe@int_colData@listData$spatialData$col)
  }

  spaObj <- createSpaTalk(st_data = st_data,
                          st_meta = st_meta,
                          species = species,
                          if_st_is_sc = is_cell,
                          spot_max_cell = 30)
  spaObj <- dec_celltype(spaObj, sc_data, sc_celltype)
  spaObj <- dec_cci_all(spaObj)
  return(spaObj)
}
