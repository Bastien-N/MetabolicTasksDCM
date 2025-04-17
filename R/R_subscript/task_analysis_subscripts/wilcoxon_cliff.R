wilcoxon_cliff <- function(dat,groups,treat,alpha = 0.05,correction = "fdr"){
  test <- require("effsize")
  test <- require("purrr") & test
  if(!test){
    stop("!! Not all packages could be loaded")
  }
  if (length(unique(groups)) != 2){
    stop("!! Wrong number of groups")
  }
  
  control <- as.data.frame(dat[groups != treat,])
  treatment <- as.data.frame(dat[groups == treat,])
  results <- map2(control,treatment,function(x,y){
      wilcox <- wilcox.test(x,y,conf.level = 1-alpha,conf.int = T)
      cd <- cliff.delta(y,x,conf.level = 1-alpha)
      res <- c(wilcox$p.value,cd$estimate,log2(mean(y)/mean(x)))
      names(res) <- c("p.value","cliff.delta","logFC")
      return(res)
    })
  results <- as.data.frame(reduce(results,rbind))
  rownames(results) <- colnames(control)
  
  results$q.value <- p.adjust(results$p.value,method = correction)
  return(results)
}
