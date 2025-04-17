#-------------------------------------------#
#   4- Comparing DCM clusters with donors   #
#-------------------------------------------#
#Last run on 12-09-2023
##- (1) Setup -#
set.seed(20230525)
#-- LocalDir definition
localDir <- "../Results/4_MAGNET_clusts_v_NF/"
if(!dir.exists(localDir)){ dir.create(localDir) }
setwd(localDir)
if(!dir.exists("Tables")){ dir.create("Tables") }
if(!dir.exists("Plots")){ dir.create("Plots")}
#-- Libraries -#
library(openxlsx)
library(ggplot2)
library(tidyr)
library(purrr)
library("stringr")
library(svglite)
library(pcaMethods)
library(glue)
library(ggprism)


source("../../R/R_subscript/task_analysis_load.R")




#-- Data -#
celScores <- t(read.table("../2_MAGNET_NFvDCM/Cellfie_scores/cellfie_score.txt",sep = ",",header = FALSE))
celScores_binary <- t(read.table("../2_MAGNET_NFvDCM/Cellfie_scores/cellfie_score_binary.txt",sep = ",",header = FALSE))
phenoData <- read.table("../0_Clean_Data/magnet_phenoData_complete.csv",sep = ",",header = T)
taskReport <- read.table("../0_Clean_Data/task_report.txt",
                         sep = ",",header = F)
taskReport <- taskReport[!taskReport$V1 %in% c("277","258"),]

TOI <- read.xlsx("../../Data/MetabolicTasks_closer look choice.xlsx",sheet = 1)
TOI <- TOI[!TOI$ID %in% c("277","258"),]
TOI[TOI$ID == "252",]$Task.of.interest <- TRUE 

#-- Data reformatting -#
phenoData$Hypertension[is.na(phenoData$Hypertension)] <- "No"
phenoData$Diabetes[is.na(phenoData$Diabetes)] <- "No"

colnames(taskReport) <- c("Task_ID","Subsystem","Source","Task_name","Failed")
taskReport <- taskReport[1:(nrow(taskReport)-2),]
taskReport$Task_ID <- paste0("Task_",taskReport$Task_ID)

taskReport$Failed <- sapply(taskReport$Failed,function(x){
  switch(x,
         "true" = F,
         T)
})

colnames(celScores) <- colnames(celScores_binary) <- taskReport$Task_ID
rownames(celScores) <- rownames(celScores_binary) <- phenoData$sampleID

minusOneScore <- celScores_binary[1,] == -1
cat(paste0(sum(taskReport$Failed)," tasks do not have score due to having failed, \n",
           sum(minusOneScore) - sum(taskReport$Failed)," tasks do not have score for other reasons, ",
           "such as not having essential rxns with associated GPR"))

celScores <- celScores[,!minusOneScore]
celScores_binary <- celScores_binary[,!minusOneScore]
taskReport <- taskReport[!minusOneScore,]
TOI <- TOI[!minusOneScore,]

mixedScore <- celScores * celScores_binary
onePct <- round(0.01*nrow(celScores_binary))
ubiquitous <- colSums(celScores_binary) <= onePct| colSums(celScores_binary) == nrow(celScores_binary) 

mixedScore_reduced <- mixedScore[,!ubiquitous]
celScores_reduced <- celScores[,!ubiquitous]
taskReport_reduced <- taskReport[!ubiquitous,]


clustTable <- read.table("../3_MAGNET_DCM_clustering/Tables/DCM_clusts_from_DCM.txt",sep = "\t",header = T)

phenoData$Cluster <- rep("Control",nrow(phenoData))
for(i in 1:nrow(clustTable)){
  phenoData$Cluster[phenoData$ID == clustTable$Sample[i]] <- clustTable$Cluster[i]
}

##- (2) PCA
pcaRes <- pcaMethods::pca(mixedScore_reduced,nPcs = 10)

plotData <- data.frame(PC1 = pcaRes@scores[,1],
                       PC2 = pcaRes@scores[,2],
                       PC3 = pcaRes@scores[,3])
plotData <- cbind(plotData,phenoData)
colorbygroup <- c("Control" = "#D55E00","Cluster_1" = "#0072B2","Cluster_2" = "#009E73")

ggplot(plotData,aes(PC1,PC2,color = Cluster))+
  geom_point()+
  scale_color_manual(values = colorbygroup)+
  theme_cardio()+
  theme(plot.background = element_rect(fill = "#ECFCFF",color = "#ECFCFF"),
        legend.background = element_rect(fill = "#ECFCFF"))+
  ylab(glue("PC2 ({round(pcaRes@R2[2]*100,1)}%)"))+
  xlab(glue("PC1 ({round(pcaRes@R2[1]*100,1)}%)"))

ggsave("Plots/PCA_point.svg",width = 12,height = 10,units = "cm")

##- (3) Hypothesis testing
#-- PhenoData
colnames(phenoData)[4] <- "ethnicity"
test_categorical_clusts <- chisq_clean(phenoData$Cluster[phenoData$Cluster %in% c("Cluster_1","Cluster_2")],
                                       phenoData[phenoData$Cluster %in% c("Cluster_1","Cluster_2"),c(3,4,8,9)])
test_categorical_clust1 <- chisq_clean(phenoData$Cluster[phenoData$Cluster %in% c("Cluster_1","Control")],
                                       phenoData[phenoData$Cluster %in% c("Cluster_1","Control"),c(3,4,8,9)])
test_categorical_clust2 <- chisq_clean(phenoData$Cluster[phenoData$Cluster %in% c("Cluster_2","Control")],
                                       phenoData[phenoData$Cluster %in% c("Cluster_2","Control"),c(3,4,8,9)])

test_categorical <- rbind(test_categorical_clust1,test_categorical_clust2)
test_categorical$Cluster <- c(rep("Cluster 1",4),rep("Cluster 2",4))
test_categorical$q.value <- p.adjust(test_categorical$p.value,method = "BH")
write.table(test_categorical,"Tables/chisq_pheno_clusters_v_control.txt",
            sep = "\t",row.names = F,col.names = T)
test_continuous_clusts <- wilcoxon_cliff(phenoData[phenoData$Cluster %in% c("Cluster_1","Cluster_2"),c(5,10,11)],
                                         groups = phenoData$Cluster[phenoData$Cluster %in% c("Cluster_1","Cluster_2")],
                                         "Cluster_1")
test_continuous_clust1 <- wilcoxon_cliff(phenoData[phenoData$Cluster %in% c("Cluster_1","Control"),c(5,10,11)],
                                   groups = phenoData$Cluster[phenoData$Cluster %in% c("Cluster_1","Control")],
                                   "Cluster_1")
test_continuous_clust2 <- wilcoxon_cliff(phenoData[phenoData$Cluster %in% c("Cluster_2","Control"),c(5,10,11)],
                                   groups = phenoData$Cluster[phenoData$Cluster %in% c("Cluster_2","Control")],
                                   "Cluster_2")

test_continuous <- rbind(test_continuous_clust1,test_continuous_clust2)
test_continuous$q.value <- p.adjust(test_continuous$p.value,method = "BH")
write.table(test_continuous,"Tables/mann_whitney_pheno_clusters_v_control.txt",
            sep = "\t",row.names = F,col.names = T)

#-- Tasks
# tmpBinary <- celScores_binary
# tmpBinary[tmpBinary == 0] <- 1e-05
# mixedScore2 <- celScores * tmpBinary

test_tasks_clust1 <- wilcoxon_cliff(celScores[phenoData$Cluster != "Cluster_2",],
                                 phenoData$Cluster[phenoData$Cluster != "Cluster_2"],
                                 "Cluster_1")
test_tasks_clust1$Comp <- c("Control_vs_Cluster1")
test_tasks_clust1$Task_ID <- row.names(test_tasks_clust1)
test_tasks_clust1$Task_name <- taskReport$Task_name

test_tasks_clust1$TOI <- TOI$Task.of.interest
test_tasks_clust1$In_Clustering <- !ubiquitous

test_tasks_clust2 <- wilcoxon_cliff(celScores[phenoData$Cluster != "Cluster_1",],
                                    phenoData$Cluster[phenoData$Cluster != "Cluster_1"],
                                    "Cluster_2")
test_tasks_clust2$Comp <- c("Control_vs_Cluster2")
test_tasks_clust2$Task_ID <- row.names(test_tasks_clust2)
test_tasks_clust2$Task_name <- taskReport$Task_name


test_tasks_clust2$TOI <- TOI$Task.of.interest
test_tasks_clust2$In_Clustering <- !ubiquitous

test_tasks_clust <- rbind(test_tasks_clust1,test_tasks_clust2)
write.table(test_tasks_clust,"Tables/mann_whitney_tasks_clusters_v_control.txt",
            sep = "\t",row.names = F,col.names = T)


# #- What if continuous
# 
# test_tasks_clust1_cont <- wilcoxon_cliff(celScores[phenoData$Cluster != "Cluster_2",],
#                                     phenoData$Cluster[phenoData$Cluster != "Cluster_2"],
#                                     "Cluster_1")
# test_tasks_clust1_cont$Comp <- c("Control_vs_Cluster1")
# test_tasks_clust1_cont$Task_ID <- row.names(test_tasks_clust1)
# test_tasks_clust1_cont$Task_name <- taskReportReduced$Task_name
# 
# test_tasks_clust2_cont <- wilcoxon_cliff(celScores[phenoData$Cluster != "Cluster_1",],
#                                     phenoData$Cluster[phenoData$Cluster != "Cluster_1"],
#                                     "Cluster_2")
# test_tasks_clust2_cont$Comp <- c("Control_vs_Cluster2")
# test_tasks_clust2_cont$Task_ID <- row.names(test_tasks_clust2)
# test_tasks_clust2_cont$Task_name <- taskReportReduced$Task_name
# 
# test_tasks_clust_cont <- rbind(test_tasks_clust1_cont,test_tasks_clust2_cont)
# 
# table(test_tasks_clust$q.value < 0.05 , test_tasks_clust_cont$q.value < 0.05)
##- (4) Vizualization
#-- PhenoData
#categorical
for(i in 1:nrow(test_categorical_clust1)){
  if(test_categorical$q.value[i] < 0.05 | test_categorical$q.value[i+nrow(test_categorical_clust1)] < 0.05){
    
    plotData <- data.frame(Var = phenoData[,test_categorical_clust1$Variable[i]],
                           Cluster = phenoData$Cluster)
    plotpval <- data.frame(group1 = c("Cluster_1","Cluster_2"),
                           group2 = c("Control","Control"),
                           label = c(glue("FDR = {format(test_categorical$q.value[i],digits = 3,scientific = -2)}"),
                                     glue("FDR = {format(test_categorical$q.value[i+nrow(test_categorical_clust1)],digits = 3,scientific = -2)}")),
                           y.position = c(1.05,1.01))
    fontToUse <- c("plain","plain")
    fontToUse[c(test_categorical$q.value[i],test_categorical$q.value[i+nrow(test_categorical_clust1)]) < 0.05] <- "bold"
    
    lineToUse <- c("dashed","dashed")
    lineToUse[c(test_categorical$q.value[i],test_categorical$q.value[i+nrow(test_categorical_clust1)]) < 0.05] <- "solid"
    ggplot(plotData,aes(x = Cluster))+
      geom_bar(stat = "count",position = "fill",aes(fill = Var))+
      theme_cardio()+
      add_pvalue(plotpval[1,],tip.length = 0.1,fontface = fontToUse[1],linetype = lineToUse[1])+
      add_pvalue(plotpval[2,],tip.length = 0.1,fontface = fontToUse[2],linetype = lineToUse[2])+
      ylab("Proportion")+
      xlab("Group")+
      labs(fill = test_categorical_clust1$Variable[i])
    ggsave(paste0("Plots/clust_v_control_",test_categorical_clust1$Variable[i],".svg"),
           height = 10,width = 12,units = "cm")
  }
}
#coninuous
for(i in 1:nrow(test_continuous_clust1)){
  if(test_continuous$q.value[i] < 0.05 | test_continuous$q.value[i+nrow(test_continuous_clust1)] < 0.05){
    
    plotData <- data.frame(Var = phenoData[,rownames(test_continuous_clust1)[i]],
                           Cluster = phenoData$Cluster)
    plotData <- plotData[!is.na(plotData[,1]),]
    plotpval <- data.frame(group1 = c("Cluster_1","Cluster_2"),
                           group2 = c("Control","Control"),
                           label = c(glue("FDR = {format(test_continuous$q.value[i],digits = 3,scientific = -2)}"),
                                     glue("FDR = {format(test_continuous$q.value[i+nrow(test_continuous_clust1)],digits = 3,scientific = -2)}")),
                           y.position = c(c(1.05,1.01)*max(plotData$Var,na.rm = T)))
    fontToUse <- c("plain","plain")
    fontToUse[c(test_categorical$q.value[i],test_categorical$q.value[i+nrow(test_categorical_clust1)]) < 0.05] <- "bold"
    
    lineToUse <- c("dashed","dashed")
    lineToUse[c(test_categorical$q.value[i],test_categorical$q.value[i+nrow(test_categorical_clust1)]) < 0.05] <- "solid"
    ggplot(plotData,aes(Cluster,Var))+
      geom_boxplot(aes(fill = Cluster))+
      geom_violin(fill = NA,linewidth = 0.1)+
      theme_cardio()+
      add_pvalue(plotpval[1,],tip.length = 0.008,fontface = fontToUse[1],linetype = lineToUse[1])+
      add_pvalue(plotpval[2,],tip.length = 0.008,fontface = fontToUse[2],linetype = lineToUse[2])+
      ylab(rownames(test_continuous_clust1)[i])+
      xlab("Group")
    ggsave(paste0("Plots/clust_v_control_",rownames(test_continuous_clust1)[i],".svg"),
           height = 10,width = 12,units = "cm")
  }
}

#-- Tasks

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
for(i in 1:nrow(test_tasks_clust1)){
  if(test_tasks_clust1$q.value[i] < 0.05 | test_tasks_clust2$q.value[i] < 0.05 ){
    if((test_tasks_clust1$q.value[i] < 0.05 & test_tasks_clust2$q.value[i] < 0.05)){
      if(sign(test_tasks_clust1$cliff.delta[i]) == sign(test_tasks_clust1$cliff.delta[i])){
        direction <- "similar"
      }else{
        direction <- "opposite"
      }
    }else{
      direction <- "different"
    }
    
    plotData <- data.frame(CellFie_score = celScores[,i],
                           Group = factor(phenoData$Cluster,levels = c("Control","Cluster_1","Cluster_2")))
    plotSigData <- data.frame(group1 = c("Control","Control"),
                              group2 = c("Cluster_1","Cluster_2"),
                              label = c(glue("FDR = {format(test_tasks_clust1$q.value[i],digits = 3,scientific = -2)}"),
                                        glue("FDR = {format(test_tasks_clust2$q.value[i],digits = 3,scientific = -2)}")),
                              y.position = c(1.01,1.08)*max(plotData[,1]))
    
    fontToUse <- c("plain","plain")
    fontToUse[c(test_tasks_clust1$q.value[i],test_tasks_clust2$q.value[i]) < 0.05] <- "bold"
    
    lineToUse <- c("dashed","dashed")
    lineToUse[c(test_tasks_clust1$q.value[i],test_tasks_clust2$q.value[i]) < 0.05] <- "solid"
    
    ggplot(plotData,aes(Group,CellFie_score))+
      geom_boxplot(aes(fill = Group))+
      geom_violin(fill = NA,linewidth = 0.1)+
      theme_cardio()+
      add_pvalue(plotSigData[1,],tip.length = 0.007,fontface = fontToUse[1],linetype = lineToUse[1])+
      add_pvalue(plotSigData[2,],tip.length = 0.007,fontface = fontToUse[2],linetype = lineToUse[2])+
      ylab("CellFie score")+
      ggtitle(paste0(rownames(test_tasks_clust1)[i],': ',formatedTaskName[i]))
    ggsave(paste0("Plots/clust_v_control_tasks_",direction,"_",test_tasks_clust1$Task_ID[i],".svg"),
           height = 10,width = 12,units = "cm")
  }
  
}

save.image("../../Data/saved_wkspaces/4_end.RData")

# 
# ################
# i=38
# 
# plotData <- data.frame(CellFie_score = mixedScore2[,i],
#                        Group = factor(phenoData$Cluster,levels = c("Control","Cluster_1","Cluster_2")))
# plotSigData <- data.frame(group1 = c("Control","Control"),
#                           group2 = c("Cluster_1","Cluster_2"),
#                           label = c(glue("FDR = {format(test_tasks_clust1$q.value[i],digits = 3,scientific = -2)}"),
#                                     glue("FDR = {format(test_tasks_clust2$q.value[i],digits = 3,scientific = -2)}")),
#                           y.position = c(1.01,1.08)*max(plotData[,1]))
# 
# fontToUse <- c("plain","plain")
# fontToUse[c(test_tasks_clust1$q.value[i],test_tasks_clust2$q.value[i]) < 0.05] <- "bold"
# 
# lineToUse <- c("dashed","dashed")
# lineToUse[c(test_tasks_clust1$q.value[i],test_tasks_clust2$q.value[i]) < 0.05] <- "solid"
# 
# colorbygroup <- c("Control" = "#D55E00","Cluster_1" = "#0072B2","Cluster_2" = "#009E73")
# 
# ggplot(plotData,aes(Group,CellFie_score))+
#   geom_boxplot(aes(fill = Group))+
#   scale_fill_manual(values = colorbygroup)+
#   geom_violin(fill = NA,linewidth = 0.1)+
#   theme_cardio()+
#   theme(plot.background = element_rect(fill = "#ECFCFF",color = "#ECFCFF"),
#         legend.background = element_rect(fill = "#ECFCFF"))+
#   add_pvalue(plotSigData[1,],tip.length = 0.007,fontface = fontToUse[1],linetype = lineToUse[1])+
#   add_pvalue(plotSigData[2,],tip.length = 0.007,fontface = fontToUse[2],linetype = lineToUse[2])+
#   ylab("CellFie score")+
#   ggtitle(paste0(rownames(test_tasks_clust1)[i],': ',formatedTaskName[i]))
# 
# ggsave(paste0("../7_Task_gene_examination/clust_v_control_tasks_",test_tasks_clust1$Task_ID[i],".svg"),
#        height = 12,width = 10,units = "cm")
