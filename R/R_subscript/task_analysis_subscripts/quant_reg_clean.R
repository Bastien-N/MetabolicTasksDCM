quant_reg_clean <- function(scores,dat){
  require("quantreg")
  require("purrr")
  tasks <- colnames(scores)
  
  qRes <- sapply(tasks,function(ta){
    dat2 <- dat
    dat2$task <- scores[tasks == ta]
    colnames(dat2)[ncol(dat2)] <- ta
    
    linmodel <- quantreg::rq(as.formula(paste0(ta," ~ .")),data = dat2)
    modelsummary <- summary(linmodel, se = "nid")
    res <- data.frame(Task = ta,Variable = rownames(modelsummary[["coefficients"]])[-1],
                      modelsummary[["coefficients"]][-1,c(1,4)])
    colnames(res)[4] <- "p.value"
    return(res)
  },simplify = F)
  qRes <- reduce(qRes,rbind)
  row.names(qRes) <- NULL
  qRes$q.value <- p.adjust(qRes[,4],"fdr")
  return(qRes)
  
  
}
