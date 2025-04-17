chisq_clean <- function(clusts,phenoDat,correct = TRUE,p = rep(1/length(x), length(x)), rescale.p = FALSE,
                        simulate.p.value = FALSE, B = 2000){
  
  pval <- rep(0,ncol(phenoDat))
  OR <- rep(0,ncol(phenoDat))
  
  for(i in 1:ncol(phenoDat)){
    if(length(unique(unlist(phenoDat[,i]))) == 1){
      pval[i] <- NA
      OR[i] <- NA
    }else{
      chires <- chisq.test(clusts,phenoDat[,i])
      pval[i] <- chires$p.value
      OR[i] <- (chires$observed[3]/chires$observed[1])/(chires$observed[4]/chires$observed[2])
    }
    
  }
  res <- data.frame(Variable = colnames(phenoDat),
                    p.value = pval,
                    odds.ratio = OR)
  return(res)
  
}

