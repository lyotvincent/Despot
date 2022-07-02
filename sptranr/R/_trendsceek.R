library(trendsceek)
# trendsceek software has some bugs in `post2pp` function, for the update of
# R Package: sf, spastat. And its latest update was 4 years ago.So we use it
# in `Seurat` as FindSpatiallyVariables(method = "markvariogram"). They are the same
# in theory.

library(Seurat)
library(cowplot)

Preprocess_Seurat <- function(seu){
  seu <- SCTransform(seu, assay = "Spatial",
                         return.only.var.genes = FALSE, verbose = FALSE)
  return(seu)
}

Estimation_trendsceek <- function(seu){
  seu <- FindSpatiallyVariableFeatures(seu, assay = "SCT", features = VariableFeatures(seu)[1:1000],
                                         selection.method = "moransi")
  top.features <- head(SpatiallyVariableFeatures(seu, selection.method = "markvariogram"), 6)
  loc <- data.frame(row = seu@images$slice1@coordinates$row,
                    col = seu@images$slice1@coordinates$col,
                    row.names = colnames(seu))
  data <- seu@assays$SCT@data
  fv <- RunMarkVario(loc, data0)
  SpatialFeaturePlot(seu, features = top.features, ncol = 3, alpha = c(0.1, 1))
}
