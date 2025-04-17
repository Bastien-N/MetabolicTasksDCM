fpc_diagnostic_parallel <- function(dat,B,nClust,seed = F){
  require(parallel)
  require(purrr)
  cl <- makeCluster(detectCores()-1)
  clusterExport(cl,c("dat","B","nClust","seed"),envir = environment())
  datFPC <- parSapply(cl,nClust,function(n){
    if (seed != F){set.seed(seed)}
    fpc::clusterboot(dat,B = B,distances = T,clustermethod = fpc::disthclustCBI,method = "ward.D2",scaling = FALSE,
                     k = n,
                     bootmethod = "boot")
    
  },simplify = FALSE)  
  bootMeandf <- sapply(datFPC,function(x){
    x$bootmean
  })
  bootMeandf <- reduce(bootMeandf,`c`)
  
  plotData <-  data.frame(Mean_jaccard = c(bootMeandf),
                          nClust = c(rep(nClust,nClust)),
                          Method = c(rep("Boot",length(bootMeandf))))
  stopCluster(cl)
  return(list(FPC_res = datFPC,
              plotData = plotData))
}