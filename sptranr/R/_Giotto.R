source("sptranr/R/_Loading.R")
Check_Load_GithubPackages("Giotto", "RubD/Giotto")

Preprocess_Giotto <- function(gobj){
  Giotto::filterCombinations(gobj,
                     expression_thresholds = c(1),
                     gene_det_in_min_cells = c(20, 20, 50, 50),
                     min_det_genes_per_cell = c(0, 200, 0, 200))
  # filter and normalize
  gobj <- Giotto::filterGiotto(gobject = gobj,
                              expression_threshold = 1,
                              gene_det_in_min_cells = 50,
                              min_det_genes_per_cell = 0,
                              expression_values = c('raw'),
                              verbose = T)
  gobj <- Giotto::normalizeGiotto(gobject = gobj, scalefactor = 6000, verbose = T)
  gobj <- Giotto::addStatistics(gobject = gobj)
  gobj <- createSpatialNetwork(gobject = gobj, minimum_k = 2, maximum_distance_delaunay = 400)
  gobj <- createSpatialNetwork(gobject = gobj, minimum_k = 2, method = 'kNN', k = 10)
  gobj <- Giotto::calculateHVG(gobject = gobj)
  return(gobj)
}

Cluster_Giotto <- function(gobj){
  # runPCA
  gobj <- Giotto::runPCA(gobject = gobj)
  gobj <- Giotto::runUMAP(gobj, dimensions_to_use = 1:10)
  gobj <- Giotto::createNearestNetwork(gobject = gobj, dimensions_to_use = 1:10, k = 15)
  gobj <- Giotto::doLeidenCluster(gobject = gobj, resolution = 0.4, n_iterations = 1000)
  return(gobj)
}

Estimate_Giotto <- function(gobj){
  spatialgenes <- binSpect(gobj)
  return(spatialgenes)
}

Deconv_Giotto <- function(gobj, sce, spatialgenes){

  ext_spatial_genes <- spatialgenes[1:100]$genes
  # 检测空间相关的基因
  spat_cor_netw_DT <- detectSpatialCorGenes(gobj,
                                           method = 'network', spatial_network_name = 'Delaunay_network',
                                           subset_genes = ext_spatial_genes)
  # 2. cluster correlation scores
  spat_cor_netw_DT <- clusterSpatialCorGenes(spat_cor_netw_DT,
                                             name = 'spat_netw_clus',
                                             k = 8)


  cluster_genes_DT <- showSpatialCorGenes(spat_cor_netw_DT,
                                          use_clus_name = 'spat_netw_clus',
                                          show_top_genes = 1)
  cluster_genes <- cluster_genes_DT$clus
  names(cluster_genes) <- cluster_genes_DT$gene_ID
  cluster_genes
  gobj <- createMetagenes(gobj, gene_clusters = cluster_genes, name = 'cluster_metagene')
  hmrf_folder = paste0(temp_dir,'11_HMRF/')
  HMRF_spatial_genes <- doHMRF(gobject = gobj,
                              expression_values = 'scaled',
                              spatial_genes = ext_spatial_genes,
                              spatial_network_name = 'Delaunay_network',
                              k = 8,
                              betas = c(28,2,2),
                              output_folder = paste0(hmrf_folder, '/', 'Spatial_genes_brain/SG_top100_k8_scaled'))
  myviewHMRFresults2D <- edit(viewHMRFresults2D)
  environment(myviewHMRFresults2D) <- environment(viewHMRFresults2D)
  debug(myviewHMRFresults2D)
  ?viewHMRFresults2D
  viewHMRFresults2D(gobject = gobj,
                    HMRFoutput = HMRF_spatial_genes,
                    k = 8, betas_to_view = 2,
                    point_size = 2)
  for(i in seq(28, 30, by = 2)) {
    myviewHMRFresults2D(gobject = gobj,
                        HMRFoutput = HMRF_spatial_genes,
                        k = 8, betas_to_view = i,
                        point_size = 2)
  }

  myaddHMRF <- edit(addHMRF)
  environment(myaddHMRF) <- environment(addHMRF)
  gobj <- myaddHMRF(gobject = gobj,
                          HMRFoutput = HMRF_spatial_genes,
                          k = 8, betas_to_add = c(28),
                          hmrf_name = 'HMRF')
  giotto_colors = getDistinctColors(8)
  spatPlot(gobject = gobj, cell_color = 'HMRF_k8_b.28',
           point_size = 3, coord_fix_ratio = 1, cell_color_code = giotto_colors)
}
