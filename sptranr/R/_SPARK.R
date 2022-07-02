# cost a lot of memory and time, not recommend.
Estimation_SPARK <- function(spark){
  ## Estimating Parameter Under Null
  spark <- spark.vc(spark,
                    covariates = NULL,
                    lib_size = spark@lib_size,
                    num_core = 5,
                    verbose = F)

  ## Calculating pval
  spark <- spark.test(spark,
                      check_positive = T,
                      verbose = F)
  return(spark)
}


Estimation_SPARKX <- function(sce){
  library(SPARK)
  counts <-sce@assays@data@listData$counts
  rownames(counts) <- rownames(sce)
  colnames(counts) <- colnames(sce)
  location <- data.frame(x = sce@colData@listData$row,
                         y = sce@colData@listData$col)
  rownames(location) <- colnames(counts)
  sparkX <- sparkx(counts,location,
                   numCores=1,option="mixture")
  return(sparkX)
}



