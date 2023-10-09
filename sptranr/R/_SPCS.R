source("sptranr/R/_Loading.R")
Check_Load_BiocPackages("Matrix")
Check_Load_BiocPackages("rsvd")
Check_Load_BiocPackages("factoextra")
Check_Load_BiocPackages("foreach")
Check_Load_BiocPackages("dplyr")
Check_Load_BiocPackages("doParallel")

#Pattern neighbors

calcRandPCA<-function(tmat, d=10, seed=42)
{
  # @param tmat A non-negative matrix with samples by features
  # @return A matrix with features by samples projected on PCA space
  set.seed(seed)

  rpca_obj <- rpca(tmat, k=d, center=T, scale=F, retx=T, p=10, q=7)
  rpca_obj$x
}

calcPDist <- function(tmat)
{
  get_dist(tmat,method = "pearson")
}

calcPContribution<-function(distMat)
{
  distMat<-as.matrix(distMat)
  rmat<-exp(-4.5*(distMat)^2)
  rmat[distMat>1]<-0
  diag(rmat)<-0

  return(Matrix(rmat, sparse=T))
}

findPNeighbors<-function(mat.pcontribution, positions, tau.p=16, nThreads=8)
{
  if(tau.p < 0)
    stop("Parameter tau.p must be 0 or positive.")

  bars<-positions$barcodes
  ind.bars<-1:length(bars)
  names(ind.bars)<-bars

  cl=makeCluster(nThreads)
  registerDoParallel(cl)

  p.neighbors<-foreach(i=1:length(bars),.packages = c("dplyr","Matrix")) %dopar%
    {
      a<-positions[mat.pcontribution[,i]>0,]
      if(nrow(a)>0){
        b<-a %>%
          mutate(p.contribution = mat.pcontribution[ind.bars[barcodes],i])
        if(0 == tau.p || nrow(a)<tau.p)
          b %>% mutate(ps.contribution = p.contribution/(sum(p.contribution)+0.00000001))
        else
          (b %>% arrange(desc(p.contribution)))[1:tau.p,] %>% mutate(ps.contribution = p.contribution/(sum(p.contribution)+0.00000001))
      } else
        NULL
    }

  stopCluster(cl)

  names(p.neighbors)<-bars

  return(p.neighbors)

}

#Spatial neighbors

GetPointNeighbors<-function(px,py,order=2,is.keepself=F)
{

  x.range<-(px-order):(px+order)
  x.range<-x.range[x.range>=0]
  y.range<-(py-order):(py+order)
  y.range<-y.range[y.range>=0]

  x.range<-rep.int(x.range,times=length(y.range))
  x.range<-sort(x.range)

  points<-data.frame(x=x.range, y=y.range)

  points<-points %>% mutate(distance=abs(x-px)+abs(y-py)) %>%
    filter(distance <= order)

  if(!is.keepself)
    points<-points %>% filter(distance>0)

  return(points)

}

findSNeighbors<-function(positions, is.hexa=F , tau.s=2, nThreads=8)
{
  if(is.hexa)
    tau.s<-tau.s*2

  bars<-positions$barcodes
  ind.bars<-1:length(bars)
  names(ind.bars)<-bars

  cl=makeCluster(nThreads)
  registerDoParallel(cl)

  s.neighbors<-foreach(i=1:length(bars),.packages =c("dplyr","Matrix"),.export = c("GetPointNeighbors")) %dopar%
    {
      a<-inner_join(x=positions, y=GetPointNeighbors(px=positions$st.x[i],py=positions$st.y[i],order = tau.s),by=c("st.x"="x","st.y"="y"))
      if(nrow(a)>0){
        a %>%
          mutate(s.contribution=(1/distance)) %>%
          mutate(ss.contribution=s.contribution/(sum(s.contribution)+0.00000001))
      } else
        NULL
    }

  stopCluster(cl)

  names(s.neighbors)<-positions$barcodes

  return(s.neighbors)

}


#get combined smoothed expression matrix

getCSmoothedExp<-function(exp, s.neighbors, p.neighbors, alpha, beta, nThreads=8)
{
  bars<-names(p.neighbors)
  ind.bars<-1:length(p.neighbors)
  names(ind.bars)<-bars

  cl=makeCluster(nThreads)
  registerDoParallel(cl)

  s.exp<-foreach(i=1:length(bars),.combine = "cbind",.packages = c("dplyr","Matrix")) %dopar%
    {
      if(!is.null(s.neighbors[[i]]))
        exp[,ind.bars[s.neighbors[[i]]$barcodes]] %*% Matrix(s.neighbors[[i]]$ss.contribution,sparse=T)
      else
        0
    }

  p.exp<-foreach(i=1:length(bars),.combine = "cbind",.packages = c("dplyr","Matrix")) %dopar%
    {
      if(!is.null(p.neighbors[[i]]))
        exp[,ind.bars[p.neighbors[[i]]$barcodes]] %*% Matrix(p.neighbors[[i]]$ps.contribution,sparse=T)
      else
        0
    }

  stopCluster(cl)

  exp.smoothed<-exp * (1-alpha) +
    (s.exp * beta + p.exp * (1-beta))*alpha

  return(exp.smoothed)
}

fillBlanks<-function(exp, coord, is.hexa=F, tau.s=2, filling.threshold = 0.5, nThreads = 8)
{

  x.bounds<-c(min(coord$st.x),max(coord$st.x))
  y.bounds<-c(min(coord$st.y),max(coord$st.y))
  x.availible<-x.bounds[1]:x.bounds[2]
  y.availible<-y.bounds[1]:y.bounds[2]
  x.availible<-sort(rep(x.availible,length(y.availible)))
  coord.maybe<-data.frame(st.x=x.availible,st.y=y.availible)

  if(is.hexa)
  {
    tau.s<-tau.s*2
    coord.maybe<- coord.maybe %>%
      mutate(is.exist = ((st.x %% 2) == (st.y %% 2))) %>%
      select(st.x,st.y)
  }

  coord.maybe<-suppressMessages(coord %>% right_join(coord.maybe,copy = T) %>% filter(is.na(barcodes)))
  coord.maybe$barcodes<-paste0("SUPP",1:nrow(coord.maybe))

  cl=makeCluster(nThreads)
  registerDoParallel(cl)

  s.neighbors<-foreach(i=1:nrow(coord.maybe),.packages = c("dplyr"),.export = c("GetPointNeighbors")) %dopar%
    {
      nbs.maybe<-GetPointNeighbors(px = coord.maybe$st.x[i],py = coord.maybe$st.y[i],order = tau.s)
      if(is.hexa)
        nbs.maybe<-nbs.maybe %>% filter(0 == distance %% 2)

      nbs<-coord %>% inner_join(nbs.maybe,by=c("st.x"="x","st.y"="y"))
      if(nrow(nbs)>filling.threshold * nrow(nbs.maybe))
      {
        nbs %>%
          mutate(s.contribution=(1/distance)) %>%
          mutate(ss.contribution=s.contribution/(sum(s.contribution)+0.00000001))
      } else
        NULL
    }

  names(s.neighbors)<-coord.maybe$barcodes
  s.neighbors<-s.neighbors[sapply(s.neighbors,function(x){!is.null(x)})]

  if(0==length(s.neighbors))
  {
    s.exp<-NULL
    s.neighbors<-NULL
  }
  else
  {
    s.exp<-foreach(i=1:length(s.neighbors),.combine = "cbind",.packages = c("dplyr","Matrix")) %dopar%
      {
        exp[,s.neighbors[[i]]$barcodes] %*% Matrix(s.neighbors[[i]]$ss.contribution,sparse=T)
      }
    colnames(s.exp)<-names(s.neighbors)
  }
  stopCluster(cl)

  return(list(exp=s.exp,colData=filter(coord.maybe,barcodes %in% names(s.neighbors))))
}

#Run SPCS
SPCS<-function(data.tbl, coor.tbl, is.hexa = F, is.padding =T, tau.p=16, tau.s=2, alpha=0.6, beta=0.4, filling.thres=0.5, nThreads=8)
{
  data.mat<-Matrix(as.matrix(data.tbl))

  message(paste0("SPCS Smoothing on ",ncol(data.tbl)," spots with ",nrow(data.tbl)," genes"))
  message("1) Pre-processing...",appendLF = F)

  t0<-proc.time()

  rpca.res = calcRandPCA(t(data.mat))
  pdist.mat = calcPDist(rpca.res)
  rmat = calcPContribution(pdist.mat)

  positions = coor.tbl
  colnames(positions) = c("st.x","st.y")
  positions$barcodes = rownames(positions)

  t1<-proc.time()
  message("Finished. Time Elapsed:",(t1-t0)["elapsed"],"secs (incl. sys. cost",(t1-t0)["sys.self"],"secs).")

  message("2) Detecting neighbors...\n"," --Pattern neighbors...",appendLF = F)
  t0<-proc.time()

  p.neighbors <- findPNeighbors(rmat, positions, tau.p, nThreads)

  t1<-proc.time()
  message("Finished. Time Elapsed:",(t1-t0)["elapsed"],"secs (incl. sys. cost",(t1-t0)["sys.self"],"secs).\n"," --Spatial neighbors...",appendLF = F)
  t0<-proc.time()

  s.neighbors <- findSNeighbors(positions,is.hexa, tau.s, nThreads)

  t1<-proc.time()
  message("Finished. Time Elapsed:",(t1-t0)["elapsed"],"secs (incl. sys. cost",(t1-t0)["sys.self"],"secs).")

  message("3) Calculating expression...",appendLF = F)
  t0<-proc.time()

  data.tbl.combinedsmooth = getCSmoothedExp(data.mat, s.neighbors, p.neighbors, alpha, beta, nThreads)

  t1<-proc.time()
  message("Finished. Time Elapsed:",(t1-t0)["elapsed"],"secs (incl. sys. cost",(t1-t0)["sys.self"],"secs).")

  if(is.padding)
  {
    message("4) Filling blank spots...",appendLF = F)
    t0<-proc.time()

    filled.blanks<-fillBlanks(data.tbl.combinedsmooth, positions, is.hexa, tau.s, filling.thres, nThreads)

    t1<-proc.time()
    message("Finished. Time Elapsed:",(t1-t0)["elapsed"],"secs (incl. sys. cost",(t1-t0)["sys.self"],"secs).\n",
            nrow(filled.blanks$colData),"blank spots has been filled.")
    data.tbl.combinedsmooth = as.data.frame(as.matrix(cbind(data.tbl.combinedsmooth,filled.blanks$exp)))
    attr(data.tbl.combinedsmooth,"SUPP_ColData")<-filled.blanks$colData
  }
  else
  {
    data.tbl.combinedsmooth = as.data.frame(as.matrix(data.tbl.combinedsmooth))
  }


  if (sum(is.na(data.tbl.combinedsmooth)) > 0) {
    mean.exp = sum(data.tbl.combinedsmooth[!is.na(data.tbl.combinedsmooth)]) / sum(!is.na(data.tbl.combinedsmooth))
    data.tbl.combinedsmooth[is.na(data.tbl.combinedsmooth)] = mean.exp
  }

  return(data.tbl.combinedsmooth)
}


# Functions for pre-processing of the ST data


# Remove samples with mostly zero genes
# Filter out samples with more than zero.cutoff*100 zero expressed genes
selectSamples <- function(data.tbl, zero.cutoff=0.99) {
  print("Filter samples")
  print(sprintf("Number of samples %d", dim(data.tbl)[2]))
  keep.idx = colSums(data.tbl==0)<=(zero.cutoff*dim(data.tbl)[1])
  data.tbl = data.tbl[,keep.idx]
  print(sprintf("Remaining number of samples %d", dim(data.tbl)[2]))
  return(list("data"=data.tbl, "idx"=keep.idx))
}


# Remove genes with mostly zeros and low variances
# Filter out genes with more than zero.cutoff*100% zero expressed samples
# Filter our genes with lowest var.cutoff*100% variances
selectGenes <- function(data.tbl, zero.cutoff=0.5, var.cutoff=0.1) {
  print("Filter genes")
  print(sprintf("Number of genes %d", dim(data.tbl)[1]))
  keep.idx1 = rowSums(data.tbl==0)<=(zero.cutoff*dim(data.tbl)[2])
  data.tbl = data.tbl[keep.idx1,]
  keep.idx = keep.idx1
  print(sprintf("Remaining number of genes %d", dim(data.tbl)[1]))
  vars = apply(data.tbl, 1, var)
  if (var.cutoff > 0) {
    keep.idx2 = vars>=quantile(vars, var.cutoff)
    data.tbl = data.tbl[keep.idx2,]
    keep.idx = keep.idx[keep.idx2]
  }
  print(sprintf("Remaining number of genes %d", dim(data.tbl)[1]))
  return(list("data"=data.tbl, "idx"=keep.idx))
}


# Normalize the library size of each sample to sum up to scale
NormalizeLibrarySize <- function(data.tbl, scale=10000) {
  return(sweep(data.tbl, 2, colSums(data.tbl), FUN="/") * scale)
}


# Decont by SPCS
Decontaminate_SPCS <- function(sce, gene.zero.cutoff = 0.7, gene.var.cutoff = 0.0,
                        tau.p=16, tau.s=2, alpha=0.6, beta=0.4, is.hexa=T,
                        nThreads=7){
  data <- sce@assays@data@listData$counts
  rownames(data) <- rownames(sce)
  colnames(data) <- colnames(sce)
  coord <- data.frame(row = sce@colData@listData$row,
                      col = sce@colData@listData$col,
                      row.names = sce@colData@listData$spot)
  filtered.genes <- selectGenes(data, gene.zero.cutoff, gene.var.cutoff)

  data <- filtered.genes$data
  coord <- coord[colnames(data),]

  data <- NormalizeLibrarySize(data)
  data <- log2(data+1)

  data.mat<-Matrix(as.matrix(data))

  rpca.res <- calcRandPCA(t(data.mat))
  pdist.mat <- calcPDist(rpca.res)
  rmat <- calcPContribution(pdist.mat)

  positions <- coord
  colnames(positions) <- c("st.x","st.y")
  positions$barcodes <- rownames(positions)

  p.neighbors <- findPNeighbors(rmat, positions, tau.p, nThreads)
  s.neighbors <- findSNeighbors(positions, tau.s, nThreads)

  data.combinedsmooth <- getCSmoothedExp(data.mat, s.neighbors, p.neighbors, alpha, beta, nThreads)
  data.combinedsmooth <- as.data.frame(as.matrix(data.combinedsmooth))

  # data.combinedsmooth <- SPCS(data, coord, tau.p, tau.s, alpha, beta, nThreads)

  SPCS_obj <- 2^data.combinedsmooth - 1
  SPCS_obj <-as(as.matrix(SPCS_obj), "dgCMatrix")
  return(SPCS_obj)
}

Save_smd_from_SPCS <- function(smdFile, SPCS_obj){
  h5data <- "SPCS_mat"

  rhdf5::h5createGroup(smdFile, h5data)
  # create barcodes
  Create_smd_array1d(smdFile,
                     arr = colnames(SPCS_obj),
                     smdloc = paste0(h5data, "/barcodes"),
                     mode = "character")

  # create data, indices, indptr, shape
  Create_smd_array1d(smdFile,
                     arr = SPCS_obj@x,
                     smdloc = paste0(h5data, "/data"),
                     mode = "double")
  Create_smd_array1d(smdFile,
                     arr = SPCS_obj@i,
                     smdloc = paste0(h5data, "/indices"),
                     mode = "integer")
  Create_smd_array1d(smdFile,
                     arr = SPCS_obj@p,
                     smdloc = paste0(h5data, "/indptr"),
                     mode = "integer")
  Create_smd_array1d(smdFile,
                     arr = SPCS_obj@Dim,
                     smdloc = paste0(h5data, "/shape"),
                     mode = "integer")

  # create features
  rhdf5::h5createGroup(smdFile, paste0(h5data, "/features"))
  Create_smd_array1d(smdFile,
                     arr = rownames(SPCS_obj),
                     smdloc = paste0(h5data, "/features/name"),
                     mode = "character")
  rhdf5::h5createGroup(smdFile, paste0(h5data, "/features/is_HVG"))

  # create `idents`
  rhdf5::h5createGroup(smdFile, paste0(h5data, "/idents"))
  # create `deconv`
  rhdf5::h5createGroup(smdFile, paste0(h5data, "/deconv"))
  # create `benchmark`
  rhdf5::h5createGroup(smdFile, paste0(h5data, "/benchmark"))
}

