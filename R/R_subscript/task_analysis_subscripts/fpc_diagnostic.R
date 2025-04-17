fpc_diagnostic <- function(dat,B,nClust,seed = F,linkage = "ward.D2"){
  require(purrr)
  datFPC <- sapply(nClust,function(n){
    if (seed != F){set.seed(seed)}
    fpc::clusterboot(dat,B = B,distances = T,clustermethod = fpc::hclustCBI,method = linkage,scaling = FALSE,
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
  return(list(FPC_res = datFPC,
              plotData = plotData))
}