recursive_regression <- function(x,y,alpha = 0.05){
  taskIDs <- colnames(y)
  linModels <- vector("list",length(taskIDs))
  for (i in 1:length(taskIDs)){
    lmData <- data.frame(y[,i],x)
    colnames(lmData)[1] <- taskIDs[i]
    isDone <- F
    while(is.data.frame(lmData) > 1 & !isDone){
      modelfit <- quantreg::rq(as.formula(paste0(taskIDs[i]," ~ .")),data = lmData)
      modelsummary <- summary(modelfit,se = 'nid')
      modelsummary
      sum((modelsummary$residuals)^2)
      if(sum(modelsummary$coefficients[,4] > alpha) == 0){
        res <- data.frame(row.names = 1:length(modelsummary$coefficients[-1,1]),
                          Task = taskIDs[i],
                          Variable = names(modelsummary$coefficients[-1,1]),
                          Coefficient = modelsummary$coefficients[-1,1],
                          p.Value = modelsummary$coefficients[-1,4])
        linModels[[i]] <- res
        isDone <- T
      }else if(sum(modelsummary$coefficients[,4] != 1) == 0){
        next()
      }else{
        maxpVal <- which.max(modelsummary$coefficients[-1,4])
        lmData <- lmData[,-(maxpVal+1)]
      }
    }
    
    
    
  }
  
  
}