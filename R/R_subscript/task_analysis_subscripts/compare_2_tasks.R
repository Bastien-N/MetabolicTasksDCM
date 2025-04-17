compare_2_tasks <- function(task1,task2,geneUsed_all,essentialRxns){
  geneUsed1 <- unique(unlist(geneUsed_all[essentialRxns[,task1] == 1,]))
  geneUsed2 <- unique(unlist(geneUsed_all[essentialRxns[,task2] == 1,]))
  
  Res <- list(Common_genes = geneUsed1[geneUsed1 %in% geneUsed2],
              t1bnt2 = geneUsed1[!(geneUsed1 %in% geneUsed2)],
              t2bnt1 = geneUsed2[!(geneUsed2 %in% geneUsed1)])
  names(Res)[2:3] <- paste0("Only_in_task_",c(task1,task2))
  return(Res)
}