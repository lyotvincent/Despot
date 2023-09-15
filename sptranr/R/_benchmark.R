library(aricode)
library(ROCR)
source("sptranr/R/_Seurat.R")

## defined a function to compute the precision, recall, and F value of two cluster results.
benchmark <- function(c1, c2){
  ##c1, a vector of reference labels that assigned to each cluster
  ##c2, a vector of predicted labels that assigned to each cluster
  #c1=true.label
  #c2=km.label.k
  ref.label <- as.vector(outer(c1, c1, "=="))
  pred.label <- as.vector(outer(c2, c2, "=="))
  pred <- prediction(as.numeric(pred.label),ref.label)

  # Accuracy
  acc.tmp <- performance(pred, "acc");
  acc <- as.numeric(acc.tmp@y.values[[1]][2])

  # ROC curve
  ROC.perf <- performance(pred, "tpr", "fpr");
  # ROC area under the curve
  auc.tmp <- performance(pred,"auc");
  auc <- as.numeric(auc.tmp@y.values)

  ##precision
  prec.tmp = performance(pred, "ppv")
  prec <- as.numeric(prec.tmp@y.values[[1]][2])
  ## F1-score
  f.tmp = performance(pred, "f")
  f <- as.numeric(f.tmp@y.values[[1]][2])

  faa <- list(F.score=f, AUC=auc, ACC=acc)

  ## ari, nmi, etc.
  clu <- clustComp(c1, c2)

  bk <- c(clu, faa)
  return(bk)
}

Benchmark_Clustering <- function(smdFile, h5data){
  h5_clu <- rhdf5::h5read(smdFile, paste0(h5data, '/idents'))
  methods <- attr(h5_clu, "names")
  mark <- data.frame()
  loc <- H5Fopen(smdFile)
  if(!H5Lexists(loc, 'matrix/idents/ground_truth')){
    message("no ground_truth for benchmark, stop.")
    H5close()
    return(mark)
  }
  H5close()
  if(h5data == "matrix"){
    methods <- methods[-which(methods == "ground_truth")]
    ground_truth <- h5_clu[["ground_truth"]]
  }else{
    # other cleaned data should load ground truth from raw data
    ground_truth <- rhdf5::h5read(smdFile, 'matrix/idents/ground_truth')
  }
  ground_truth[is.na(ground_truth)] <- "NA"
  ground_truth <- as.factor(ground_truth)
  for(met in methods){
    comp <- data.frame(benchmark(as.factor(h5_clu[[met]]), ground_truth),
                       row.names = met)
    mark <- rbind(mark, comp)
  }
  return(mark)
}

Save_smd_from_Benchmark_clu <- function(smdFile, h5data, mark){
  # convert to dense Matrix
  mat_mark <- Matrix(as.matrix(mark))

  h5createGroup(smdFile, paste0(h5data, '/benchmark/cluster'))
  # save 1d weights
  Create_smd_array1d(smdFile,
                     arr = mat_mark@x,
                     smdloc = paste0(h5data, '/benchmark/cluster/value'),
                     mode = 'double')
  # save shape
  Create_smd_array1d(smdFile,
                     arr = mat_mark@Dim,
                     smdloc = paste0(h5data, '/benchmark/cluster/shape'),
                     mode = 'integer')
  # save dim names
  Create_smd_array1d(smdFile,
                     arr = rownames(mat_mark),
                     smdloc = paste0(h5data, '/benchmark/cluster/methods'),
                     mode = 'character')
  Create_smd_array1d(smdFile,
                     arr = colnames(mat_mark),
                     smdloc = paste0(h5data, '/benchmark/cluster/targets'),
                     mode = 'character')
}

moransI <- function(data, loc, features){
  filtered_data <- data[features, ]
  moran <- RunMoransI(filtered_data, loc)
  return(moran)
}

Benchmark_Estimation <- function(smdFile,imgdir, h5data){
  seu <-Load_smd_to_Seurat(smdFile, imgdir, h5data)
  # load location info
  loc <- data.frame(row = seu@images$slice1@coordinates$row,
                    col = seu@images$slice1@coordinates$col,
                    row.names = colnames(seu))
  seu <- Preprocess_Seurat(seu)
  data <- seu@assays$SCT@data

  h5_est <- rhdf5::h5read(smdFile, paste0(h5data, '/features/is_HVG'))

  methods <- attr(h5_est, "names")
  moran = list()
  for(method in methods){
    if(method == "BayesSpace"){
      # load BayesSpace HVGs
      fea <- data.frame(is_HVG = h5_est[[method]],
                              gene = rownames(seu),
                              row.names = rownames(seu))
      fea <- rownames(fea[which(fea == T),])
    }
    else if(method == "leiden"){
      # load leiden HVGs
      fea <- data.frame(is_HVG = h5_est[[method]],
                               gene = rownames(seu),
                               row.names = rownames(seu))
      fea <- rownames(fea[which(fea$is_HVG == T),])
    }
    else if(method == "SPARK"){
      # load SPARK HVGs
      fea <- data.frame(gene = h5_est[[method]]$name,
                              qval = h5_est[[method]]$qval,
                              row.names = h5_est[[method]]$name)
      fea <- rownames(fea[order(fea$qval, decreasing = FALSE)[1:500],])
    }
    else if(method == "SpaGCN"){
      # load SpaGCN HVGs
      fea <- data.frame(gene = h5_est[[method]]$name,
                               pval = h5_est[[method]]$pval)
      fea <- fea[!duplicated(fea$gene), ]$gene
    }
    else if(method == "SpatialDE"){
      # load SpatialDE HVGs
      fea <- data.frame(gene = h5_est[[method]]$name,
                                  qval = h5_est[[method]]$qval,
                                  row.names = h5_est[[method]]$name)
      fea <- rownames(fea[order(fea$qval, decreasing = FALSE)[1:500],])
    }
    else if(method == "Giotto"){
      # load SpatialDE HVGs
      fea <- data.frame(gene = h5_est[[method]]$name,
                        qval = h5_est[[method]]$qval,
                        row.names = h5_est[[method]]$name)
      fea <- rownames(fea[order(fea$qval, decreasing = FALSE)[1:500],])
    }
    fea.moran <- moransI(data, loc, fea)
    moran[[method]] <- fea.moran
  }

  return(moran)
}

Save_smd_from_Benchmark_est <- function(smdFile, h5data, moran){

  h5createGroup(smdFile, paste0(h5data, '/benchmark/estimate'))

  methods <- attr(moran, "names")

  for(method in methods){

    h5createGroup(smdFile, paste0(h5data, '/benchmark/estimate/', method))
    # save gene
    Create_smd_array1d(smdFile,
                       arr = rownames(moran[[method]]),
                       smdloc = paste0(h5data, '/benchmark/estimate/', method, "/name"),
                       mode = 'character')
    # save moransI
    Create_smd_array1d(smdFile,
                       arr = moran[[method]]$observed,
                       smdloc = paste0(h5data, '/benchmark/estimate/', method, "/moransI"),
                       mode = 'double')
    # save pvalue
    Create_smd_array1d(smdFile,
                       arr = moran[[method]]$p.value,
                       smdloc = paste0(h5data, '/benchmark/estimate/', method, "/pval"),
                       mode = 'double')
  }
}
