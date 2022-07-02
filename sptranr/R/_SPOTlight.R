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

# save spt from SPOTlight results
Save_spt_from_SPOTlight <- function(sptFile, h5data, spotlight){

  results <- as(spotlight$mat, "dgeMatrix")

  h5createGroup(sptFile, paste0(h5data, '/deconv/SPOTlight'))

  # save 1d weights
  Create_spt_array1d(sptFile,
                     arr = results@x,
                     sptloc = paste0(h5data, '/deconv/SPOTlight/weights'),
                     mode = 'double')
  # save shape
  Create_spt_array1d(sptFile,
                     arr = results@Dim,
                     sptloc = paste0(h5data, '/deconv/SPOTlight/shape'),
                     mode = 'integer')
  # save dim names
  Create_spt_array1d(sptFile,
                     arr = rownames(results),
                     sptloc = paste0(h5data, '/deconv/SPOTlight/barcodes'),
                     mode = 'character')
  Create_spt_array1d(sptFile,
                     arr = colnames(results),
                     sptloc = paste0(h5data, '/deconv/SPOTlight/cell_type'),
                     mode = 'character')
}
