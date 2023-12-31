source("sptranr/R/_Loading.R")
Check_Load_GithubPackages("spacexr", "https://github.com/dmcable/spacexr")

GenerateRef_spacexr <- function(sce){
  counts <- sce@assays@data$counts
  counts <- as(counts, "sparseMatrix")
  rownames(counts) <- rownames(sce)
  colnames(counts) <- colnames(sce)
  cell_type <- as.factor(sce@colData@listData$free_annotation)
  names(cell_type) <- colnames(sce)
  ref <- spacexr::Reference(counts,
                   cell_types = cell_type)

  return(ref)
}

Deconvolution_spacexr <- function(sr, ref){
  rctd <- create.RCTD(sr, ref,
                      max_cores = 7,
                      UMI_min = 0,
                      CELL_MIN_INSTANCE=10)
  rctd <- run.RCTD(rctd)
  return(rctd)
}

Save_smd_from_spacexr <- function(smdFile, h5data, rctd, pp_mtd="vae"){
  results <- rctd@results
  # normalize the cell type proportions to sum to 1.
  norm_weights <- normalize_weights(results$weights)
  cell_type_names <- rctd@cell_type_info$info[[2]] #list of cell type names
  if(pp_mtd == "vae"){
    mtd <- 'spacexr'
  } else if(pp_mtd == "es"){
    mtd <- 'spacexr_es'
  }
  h5createGroup(smdFile, paste0(h5data, '/deconv/', mtd))

  # save 1d weights
  Create_smd_array1d(smdFile,
                     arr = norm_weights@x,
                     smdloc = paste0(h5data, '/deconv/', mtd, '/weights'),
                     mode = 'double')
  # save shape
  Create_smd_array1d(smdFile,
                     arr = norm_weights@Dim,
                     smdloc = paste0(h5data, '/deconv/', mtd, '/shape'),
                     mode = 'integer')
  # save dim names
  Create_smd_array1d(smdFile,
                     arr = rownames(norm_weights),
                     smdloc = paste0(h5data, '/deconv/', mtd,'/barcodes'),
                     mode = 'character')
  Create_smd_array1d(smdFile,
                     arr = colnames(norm_weights),
                     smdloc = paste0(h5data, '/deconv/', mtd, '/cell_type'),
                     mode = 'character')
}

Plot_spacexr <- function(rctd, figdir){
  results <- rctd@results
  # normalize the cell type proportions to sum to 1.
  norm_weights <- normalize_weights(results$weights)
  cell_type_names <- rctd@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- rctd@spatialRNA
  resultsdir <- 'RCTD_Plots' ## you may change this to a more accessible directory on your computer.
  dir.create(resultsdir)
  # make the plots
  # Plots the confident weights for each cell type as in full_mode (saved as
  # 'results/cell_type_weights_unthreshold.pdf')
  plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
  # Plots all weights for each cell type as in full_mode. (saved as
  # 'results/cell_type_weights.pdf')
  plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights)
  # Plots the weights for each cell type as in doublet_mode. (saved as
  # 'results/cell_type_weights_doublets.pdf')
  plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet,
                       results$results_df)
  # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as
  # 'results/cell_type_occur.pdf')
  plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
  # makes a map of all cell types, (saved as
  # 'results/all_cell_types.pdf')
  plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir)
}
