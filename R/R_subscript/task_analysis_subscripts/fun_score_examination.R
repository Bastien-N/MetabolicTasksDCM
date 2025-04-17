get_rxn_score <- function(taskInd,essentialRxns,activity_scores,promiscuity){
  score <- as.matrix(activity_scores[essentialRxns[,taskInd] == 1,]) * 1/as.matrix(promiscuity[essentialRxns[,taskInd] == 1,])
  notNa <- !is.na(score[,1])
  if(sum(notNa) == 0){
    errorCondition("No essential reaction has an activity score")
  }else{
    score <- score[!is.na(score[,1]),]
    return(score)
  }
}
important_rxn_global <- function(taskInd,essentialRxns,activity_scores,promiscuity){
  scores <- get_rxn_score(taskInd,essentialRxns,activity_scores,promiscuity)
  
  avgSum <- mean(colSums(scores))
  rxnAvg <- rowMeans(scores)
  rxnSd <- apply(scores,1,sd)
  
  Res <- data.frame(Essential_rxns = row.names(scores))
  restmp <- data.frame(Rxn_avg = rxnAvg,
                       Rxn_pct = rxnAvg/avgSum,
                       Rxn_sd = rxnSd)
  Res <- cbind(Res,restmp)
  Res <- Res[order(Res[,3],decreasing = T),]
  return(Res)
}
important_rxn_by_group <- function(taskInd,groups,essentialRxns,activity_scores,promiscuity){
  scores <- get_rxn_score(taskInd,essentialRxns,activity_scores,promiscuity)
  gs <- unique(groups)
  avgSum <- rep(0,length(gs)) 
  rxnAvg <- matrix(ncol = length(gs), nrow = nrow(scores)) 
  rxnSd <- matrix(ncol = length(gs), nrow = nrow(scores)) 
  for(i in 1:length(gs)){
    avgSum[i] <- mean(colSums(scores[,groups == gs[i]]))
    rxnAvg[,i] <- rowMeans(scores[,groups == gs[i]])
    rxnSd[,i] <- apply(scores[,groups == gs[i]],1,sd)
  }
  Res <- data.frame(Essential_rxns = row.names(scores))
  for(i in 1:length(gs)){
    restmp <- data.frame(Rxn_avg = rxnAvg[,i],
                         Rxn_pct = rxnAvg[,i]/avgSum[i],
                         Rxn_sd = rxnSd[,i])
    colnames(restmp) <- paste0(colnames(restmp),"_",gs[i])
    Res <- cbind(Res,restmp)
  }
  Res <- Res[order(Res[,3],decreasing = T),]
  
  return(Res)
}

format_rxnImpTable <- function(rxnImpTable,rxnTable,geneUsed){
  require(dplyr)
  ng <- (ncol(rxnImpTable)-1)/3
  Res <- rxnImpTable
  for (n in 1:ng){
    Res[,1+(n-1)*3+2] <- paste0(round(Res[,1+(n-1)*3+2]*100,1),"%") 
  }
  genes <- rep(NA,nrow(Res))
  for(i in 1:nrow(Res)){
    genes[i] <- paste0(unique(unlist(geneUsed[rxnImpTable$Essential_rxns[i],])),collapse = "; ")
  }
  Res <- right_join(rxnTable,Res,join_by(Rxn_ID == Essential_rxns))
  Res <- Res[order(as.numeric(gsub("\\%","",Res[,6])),decreasing = T),]
  Res$Genes <- genes
  Res$ratio_DCM_NF <-Res$Rxn_avg_DCM / Res$Rxn_avg_NF
  return(Res)
}