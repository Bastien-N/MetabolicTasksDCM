#-----------------------------------#
#   3- Clustering MAGNET DCM Data   #
#-----------------------------------#
#Last run on 11-09-2023
##- (1) Setup -#
set.seed(20230525)
#-- LocalDir definition
localDir <- "../Results/3_MAGNET_DCM_clustering/"
if(!dir.exists(localDir)){ dir.create(localDir) }
setwd(localDir)
if(!dir.exists("Tables")){ dir.create("Tables") }
if(!dir.exists("Plots")){ dir.create("Plots") }
#-- Libraries -#
library(openxlsx)
library(ggplot2)
library(tidyr)
library(purrr)
library("stringr")
library(pheatmap)
library(fpc)
library(svglite)
library(pcaMethods)
library(partykit)
library(glue)
library(ggprism)
library(gtsummary)
library(dplyr)
source("../../R/R_subscript/task_analysis_load.R")


#-- Plots handling -#
createdPlots <- vector("list")

#-- Data -#
celScores <- t(read.table("Cellfie_scores/cellfie_score.txt",sep = ",",header = FALSE))
celScores_binary <- t(read.table("Cellfie_scores/cellfie_score_binary.txt",sep = ",",header = FALSE))
phenoData <- read.table("../0_Clean_Data/magnet_phenoData_DCM.csv",sep = ",",header = T)
taskReport <- read.table("../0_Clean_Data/task_report.txt",
                         sep = ",",header = F)
taskReport <- taskReport[!taskReport$V1 %in% c("277","258"),]

TOI <- read.xlsx("../../Data/MetabolicTasks_closer look choice.xlsx",sheet = 1)
TOI <- TOI[!TOI$ID %in% c("277","258"),]
TOI[TOI$ID == "252",]$Task.of.interest <- TRUE 



#-- Data reformatting -#
colnames(taskReport) <- c("Task_ID","Subsystem","Source","Task_name","Failed")
taskReport <- taskReport[1:(nrow(taskReport)-2),]
taskReport$Task_ID <- paste0("Task_",taskReport$Task_ID)

taskReport$Failed <- sapply(taskReport$Failed,function(x){
  switch(x,
         "true" = F,
         T)
})

colnames(celScores) <- colnames(celScores_binary) <- taskReport$Task_ID
rownames(celScores) <- rownames(celScores_binary) <- phenoData$ID

minusOneScore <- celScores[1,] == -1
cat(paste0(sum(taskReport$Failed)," tasks do not have score due to having failed, \n",
           sum(minusOneScore) - sum(taskReport$Failed)," tasks do not have score for other reasons, ",
           "such as not having essential rxns with associated GPR"))

celScores <- celScores[,!minusOneScore]
celScores_binary <- celScores_binary[,!minusOneScore]
taskReport <- taskReport[!minusOneScore,]
TOI <- TOI[!minusOneScore,]

mixedScore <- celScores * celScores_binary
onePct <- round(0.01*nrow(celScores_binary))
ubiquitous <- colSums(celScores_binary) <= onePct  | colSums(celScores_binary) == nrow(celScores_binary) 


mixedScore_reduced <- mixedScore[,!ubiquitous]
celScores_reduced <- celScores[,!ubiquitous]
taskReport_reduced <- taskReport[!ubiquitous,]
##- (2) Clustering -#
#-- PCA -#
pcaRes <- pcaMethods::pca(mixedScore_reduced,nPcs = ncol(mixedScore_reduced))
pcaScores <- pcaRes@scores[,pcaRes@R2 > 0.01]

#-- PCA exploratory vizualization
plotData <- data.frame(PC1 = pcaRes@scores[,1],PC2 = pcaRes@scores[,2],
                       gender = phenoData$gender,race = phenoData$race,
                       HT = phenoData$Hypertension,Diabetes = phenoData$Diabetes)

ggplot(plotData, aes(PC1,PC2,color = race))+
  geom_point()
ggsave("Plots/Exploratory/pca_race_point.png")


ggplot(plotData,aes(race,PC1))+
  geom_boxplot()
ggsave("Plots/Exploratory/pca_race_boxplot.png")

  
#-- FPC run -#
B <- 10000 #Number of boostraps
nClusts <- 2:15 #Number of clusters to attempt

# fpcRes <- fpc_diagnostic_parallel(pcaScores,B,nClusts,seed = 20222023)
fpcRes <- fpc_diagnostic_parallel(dist(pcaScores),B,nClusts,seed = 20222023)
save(fpcRes,file = "fpcRes.RData")
#-- FPC Visualization -#
createdPlots$p1 <- ggplot(fpcRes$plotData,aes(x=nClust,y=Mean_jaccard))+
  geom_point()+
  geom_hline(yintercept = 0.6)+
  theme_cardio()+
  ylab("Mean JC across bootstraps")+
  xlab("Number of clusters")

plot(createdPlots$p1)
ggsave("Plots/1_FPC_plot.svg",createdPlots$p1)

#-- Clustering -#
allAbove <- rep(F,length(nClusts))
for (i in 1:length(nClusts)){
  n <- nClusts[i]
  allAbove[i] <- sum(fpcRes$plotData$nClust == n & fpcRes$plotData$Mean_jaccard > 0.6) == n
}
nClustInd <- max(which(allAbove))

cat(paste0("There are ",nClusts[nClustInd]," stable clusters"))

clusters <- fpcRes$FPC_res[[nClustInd]]$partition
clusters <- paste0("Cluster_",clusters)
clustTable <- data.frame(Sample = phenoData$ID,
                         Cluster = clusters)
#save.image("../../Data/saved_wkspaces/saved_for_abstr.RData")
load("../../Data/saved_wkspaces/saved_for_abstr.RData")

write.table(clustTable,"Tables/DCM_clusts_from_DCM.txt",sep = "\t",col.names = T,row.names = F)
clustTable <- read.table("Tables/DCM_clusts_from_DCM.txt",sep = "\t",header = T)
clusters <- clustTable$Cluster

#-- Cluster Visualization -#
# sampleAnnotation = data.frame(Tissue = phenoData$tissue,
#                               row.names = phenoData$ID)
mixedScore_reduced_scaled <- apply(mixedScore_reduced,2,function(x){x / max(x)})
svg("Plots/2_Cluster_heatmap.svg")
pheatmap(t(mixedScore_reduced_scaled),cluster_cols = fpcRes$FPC_res[[nClustInd]]$result$result,
         #annotation_col = sampleAnnotation,
         show_rownames = FALSE,
         clustering_method = "ward.D2",
         treeheight_row = 10,
         show_colnames = FALSE,method = "ward.D2",cutree_cols = nClusts[nClustInd])
dev.off(dev.cur())

png("Plots/2_Cluster_heatmap.png",width = 30,height = 25,units = "cm",res = 300)
pheatmap(t(mixedScore_reduced_scaled),cluster_cols = fpcRes$FPC_res[[nClustInd]]$result$result,
         #annotation_col = sampleAnnotation,
         show_rownames = FALSE,
         clustering_method = "ward.D2",
         treeheight_row = 10,
         show_colnames = FALSE,method = "ward.D2",cutree_cols = nClusts[nClustInd]
        )
graphics.off()
##- (2) Cluster examination -#
pcaRes <- pcaMethods::pca(mixedScore_reduced)

plotData <- data.frame(PC1 = pcaRes@scores[,1],
                       PC2 = pcaRes@scores[,2],
                       Cluster = clusters)
ggplot(plotData,aes(PC1,PC2,color = Cluster))+
  geom_point()+
  scale_color_colorblind()+
  theme_cardio()+
  ylab(glue("PC2 ({round(pcaRes@R2[2]*100,1)}%)"))+
  xlab(glue("PC1 ({round(pcaRes@R2[1]*100,1)}%)"))
ggsave("Plots/PCA_point.svg",width = 12,height = 10,units = "cm")

#-- Relationship phenoData and clusters
test_categorical <- chisq_clean(clusters,phenoData[,c(2,3,7,8)])
test_categorical$q.value 

write.table(test_categorical,"Tables/chisq_clusters_variables.txt",sep = "\t",col.names = T,row.names = F)

test_continuous <- wilcoxon_cliff(phenoData[,c(4,9,10)],groups = as.character(clusters),"Cluster_2")

write.table(test_continuous,"Tables/mann_whitney_clusters_variables.txt",sep = "\t",col.names = T,row.names = F)

phenoData |> 
  mutate(Metabotype = str_replace(clusters,"Cluster","Metabotype")) |> 
  mutate(race = str_replace(race,"AA","Afr. American")) |> 
  rename(Ethnicity = race,Age = age,Sex = gender) |> 
  gtsummary::tbl_summary(include = c(Sex,Age,Ethnicity,BMI,
                                     Diabetes,Hypertension,LVEF),
                         by = Metabotype,
                         missing = "ifany") |> 
  gtsummary::add_p() |> 
  as_gt() |> 
  gt::gtsave("Tables/Phenotype_v_clust.docx")

#-- Cluster characterization 



test_cluster <- wilcoxon_cliff(celScores,as.character(clusters),"Cluster_2",correction = "BH")
test_cluster$Task_ID <- row.names(test_cluster)
test_cluster$Task_name <- str_remove(taskReport$Task_name,"\\(.*$")
test_cluster$TOI <- TOI$Task.of.interest
test_cluster$In_Clustering <- !ubiquitous
write.table(test_cluster,"Tables/mann_whitney_cluster_tasks.txt",sep = "\t",col.names = T,row.names = F)

# 
# #Conditional forest
# cfData <- as.data.frame(mixedScore)
# cfData$Cluster <- as.factor(clusters)
# forestObj <- cforest(Cluster ~ .,data = cfData)
# 
# predTest <- predict(forestObj,cfData[,-ncol(cfData)],type = "response")
# table(predTest,cfData$Cluster)
# accuracy <- sum(predTest == cfData$Cluster)/length(predTest)
# cat(paste0("The conditionnal forest is accurate in ",round(accuracy*100,2),"% of cases"))
# 
# #Conditional variable importance
# cVarImp <- partykit::varimp(forestObj,conditional = T,nperm = 10)
# 
# taskImportance <- data.frame(Task_ID = taskReport$Task_ID[taskReport$Task_ID %in% names(cVarImp)],
#                              Task_name = taskReport$Task_name[taskReport$Task_ID %in% names(cVarImp)],
#                              Mean_acc_loss = cVarImp)
# taskImportance <- taskImportance[order(cVarImp,decreasing = T),]
# 
# ggplot(cfData,aes(Cluster,Task_274))+
#   geom_violin()

##- (3) Visualization

formatedTaskName <- taskReport$Task_name
formatedTaskName <- str_remove(formatedTaskName," \\(.*$")
for(i in 1:length(formatedTaskName)){
  if(nchar(formatedTaskName[i]) > 50){
    spaces <- str_locate_all(formatedTaskName[i]," ")
    lastSpace <- max(spaces[[1]][spaces[[1]][,1] <= 50,1])
    formatedTaskName[i] <- paste0(substr(formatedTaskName[i],1,lastSpace-1),
                                  "\n",
                                  substr(formatedTaskName[i],lastSpace+1,nchar(formatedTaskName[i])))
  }
}


for(i in 1:nrow(test_cluster)){
  if(test_cluster$q.value[i] < 0.05){
    plotData <- data.frame(score = celScores[,i],
                           Sex = phenoData$gender,
                           Age = phenoData$age,
                           Race = phenoData$race,
                           Cluster = as.character(clusters))
    
    sigTable <- data.frame(
      group1 = "Cluster_1",
      group2 = "Cluster_2",
      p.adj = test_cluster$q.value[i],
      y.position = c(max(celScores[,i])*1.05)
    )
    ggplot(plotData,aes(Cluster,score))+
      geom_boxplot(aes(fill = Cluster))+
      geom_violin(linewidth = 0.2,fill = NA)+
      #geom_point(color = "black",alpha = 0.4,position = position_jitter(height = 0,width = 0.05))+
      ylab("CellFie Score")+
      ggtitle(paste0(colnames(celScores)[i],": ",formatedTaskName[i]))+
      theme_cardio()+
      labs(fill = "Cluster")+
      xlab("Cluster")+
      add_pvalue(data = sigTable,"FDR = {format(p.adj,digits = 3,scientific = -2)}",tip.length = 0.01)
    ggsave(paste0("Plots/mann_whitney/",colnames(celScores)[i],"_clust1_v_clust2.svg"))
  }
}

ggplot(test_cluster,aes(logFC,-log(p.value))) +
  geom_point() +
  geom_hline(yintercept = -log(0.05),color = "red")
ggsave("Plots/volcano_plot_pvals.png")

ggplot(test_cluster,aes(logFC,-log(q.value))) +
  geom_point() +
  geom_hline(yintercept = -log(0.05),color = "red")
ggsave("Plots/volcano_plot_qvals.png")


save.image("../../Data/saved_wkspaces/3_end.RData")
#load("../../Data/saved_wkspaces/3_end.RData")

# #################
# probs <- seq(0,1,0.001)
# cliffs <- rep(0,length(probs))
# n <- 10000
# for(i in 1:length(cliffs)){
#   x <- c(rep(1,probs[i]*n),rep(0,(1-probs[i])*n))
#   y <- rep(0.5,n)
#   cd <- cliff.delta(x,y)
#   cliffs[i] <- cd$estimate
# }
# plot(cliffs,probs)
# 
# test <- read.table("geTMM_MAGNET_BC.csv",sep = ",",header = T)
# test2 <- read.table("magnet_phenoData_complete.csv",sep = ",",header = T)
# 
# qntGlobal <- quantile(log(1+unlist(test[,-1])),c(0.25,0.75))
# qntDCM <- quantile(log(1+unlist(test[plotData$Group == "DCM",-1])),c(0.25,0.75))
# qntControl <- quantile(log(1+unlist(test[plotData$Group == "NF",-1])),c(0.25,0.75))
# library(patchwork)
# i=2
# plotData <- data.frame(Expression = log(1+unlist(test[i,-1])),
#                        Group = test2$etiology,
#                        All = "All")
# meanGlobal <- mean(plotData$Expression)
# meanDCM <- mean(plotData[plotData$Group == "DCM",]$Expression)
# meanControl <- mean(plotData[plotData$Group == "NF",]$Expression)
# 
# g <- ggplot(plotData[plotData$Group == "DCM",],aes(Group,Expression))+
#   geom_violin()+
#   geom_hline(color = "blue",yintercept = qntDCM[1])+
#   geom_hline(color = "blue",yintercept = qntDCM[2])+
#   geom_hline(color = "red",yintercept = meanDCM)+
#   ggtitle("DCM")
# g2 <- ggplot(plotData[plotData$Group == "NF",],aes(Group,Expression))+
#   geom_violin()+
#   geom_hline(color = "blue",yintercept = qntControl[1])+
#   geom_hline(color = "blue",yintercept = qntControl[2])+
#   geom_hline(color = "red",yintercept = meanControl)+
#   ggtitle("Control")
# g3 <- ggplot(plotData,aes(All,Expression))+
#   geom_violin()+
#   geom_hline(color = "blue",yintercept = qntGlobal[1])+
#   geom_hline(color = "blue",yintercept = qntGlobal[2])+
#   geom_hline(color = "red",yintercept = meanGlobal)+
#   ggtitle("All")
# 
# g+g2+g3
# Expr <- log(1+unlist(test[i,-1]))
# thresDat <- data.frame(Group = c("All","DCM","Control"),
#                        
#                        M = c(meanGlobal,meanDCM,meanControl),
#                        p25 = c(qntGlobal[1],qntDCM[1],qntDCM[1]),
#                        p75 = c(qntGlobal[2],qntDCM[2],qntDCM[2]))
# plotData <- data.frame(Expression = c(Expr,Expr[test2$etiology == "DCM"],Expr[test2$etiology == "NF"]),
#                        X = "x",
#                        Group = c(rep("All",length(Expr)),
#                                  rep("DCM",sum(test2$etiology == "DCM")),
#                                  rep("Control",sum(test2$etiology == "NF"))))
# g4 <- ggplot(plotData,aes(X,Expression))+
#   geom_violin()+
#   
#   geom_hline(data = thresDat,aes(yintercept = p25),color = "blue")+
#   facet_wrap(vars(Group))
#   geom_hline(color = "blue",yintercept = qntGlobal[2])+
#   geom_hline(color = "red",yintercept = meanGlobal)
