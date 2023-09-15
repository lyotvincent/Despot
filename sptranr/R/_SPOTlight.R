library(SPOTlight)

# Load SCE data
Deconvolution_SPOTlight <- function(spe, sce){

  res <- SPOTlight(
    x = sce,
    y = spe,
    assay = "counts",
    groups = sce$free_annotation,
    hvg = sce@metadata$HVGs,
    mgs = sce@metadata$marker_genes,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")

  return(res)
}

# save smd from SPOTlight results
Save_smd_from_SPOTlight <- function(smdFile, h5data, spotlight, pp_mtd="raw"){

  results <- as(spotlight$mat, "dgeMatrix")

  if(pp_mtd == "raw"){
    mtd_name <- "SPOTlight"
  } else if(pp_mtd == "es"){
    mtd_name <- "SPOTlight_es"
  } else if(pp_mtd == "vae"){
    mtd_name <- "SPOTlight_vae"
  }

  h5createGroup(smdFile, paste0(h5data, '/deconv/', mtd_name))

  # save 1d weights
  Create_smd_array1d(smdFile,
                     arr = results@x,
                     smdloc = paste0(h5data, '/deconv/', mtd_name, '/weights'),
                     mode = 'double')
  # save shape
  Create_smd_array1d(smdFile,
                     arr = results@Dim,
                     smdloc = paste0(h5data, '/deconv/', mtd_name ,'/shape'),
                     mode = 'integer')
  # save dim names
  Create_smd_array1d(smdFile,
                     arr = rownames(results),
                     smdloc = paste0(h5data, '/deconv/', mtd_name ,'/barcodes'),
                     mode = 'character')
  Create_smd_array1d(smdFile,
                     arr = colnames(results),
                     smdloc = paste0(h5data, '/deconv/', mtd_name ,'/cell_type'),
                     mode = 'character')
}
