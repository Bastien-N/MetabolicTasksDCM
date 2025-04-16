#-----------------------------#
#   1- Clustering GTEx Data   #
#-----------------------------#
##- (1) Setup -#
#-- LocalDir definition
localDir <- "../Results/1_GTEx_clustering/"
setwd(localDir)
#-- Libraries -#
library(ggplot2)
library(tidyr)
library(purrr)
library("stringr")
library(pheatmap)
library(fpc)
library(svglite)
library(rstatix)
library(pcaMethods)
library(openxlsx)
library(ggpubr)
library(patchwork)
source("../../R/R_subscript/task_analysis_load.R")

#-- Plots handling -#
createdPlots <- vector("list")

#-- Data -#
celScores <- t(read.table("Cellfie_scores/cellfie_score.txt",sep = ",",header = FALSE))
celScores_binary <- t(read.table("Cellfie_scores/ cellfie_score_binary.txt",sep = ",",header = FALSE))
phenoData <- read.table("../0_Clean_Data/GTEx_phenoData.csv",sep = ",",header = T)
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
rownames(celScores) <- rownames(celScores_binary) <- phenoData$sampleID

minusOneScore <- celScores[1,] == -1
cat(paste0(sum(taskReport$Failed)," tasks do not have score due to having failed, \n",
             sum(minusOneScore) - sum(taskReport$Failed)," tasks do not have score for other reasons, ",
            "such as not having essential rxns with associated GPR"))

celScores <- celScores[,!minusOneScore]
celScores_binary <- celScores_binary[,!minusOneScore]
taskReport <- taskReport[!minusOneScore,]
TOI <- TOI[!minusOneScore,]

mixedScore <- celScores * celScores_binary

ubiquitous <- colSums(celScores_binary) == 0 | colSums(celScores_binary) == nrow(celScores_binary) 
mixedScore <- mixedScore[,!ubiquitous]

##- (2) Clustering -#
#-- PCA -#
pcaRes <- pcaMethods::pca(mixedScore,nPcs = ncol(mixedScore))
pcaScores <- pcaRes@scores[,pcaRes@R2 > 0.001]

#-- FPC run -#
B <- 10000 #Number of boostraps
nClusts <- 2:15 #Number of clusters to attempt

fpcRes <- fpc_diagnostic_parallel(dist(pcaScores),B,nClusts,20222023)
# 
# corMat <- cor(t(mixedScore))
# corDist <- as.dist(1-corMat)
# fpcTest <- fpc_diagnostic_parallel(corDist,10,nClusts,20222023)

#-- DPC Visualization -#
createdPlots$p1 <- ggplot(fpcRes$plotData,aes(x=nClust,y=Mean_jaccard))+
  geom_point()+
  geom_hline(yintercept = 0.6)+
  theme_cardio()+
  ylab("Mean JC across bootstraps")
  
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

#-- Cluster Visualization -#
sampleAnnotation = data.frame(Tissue = phenoData$tissue,
                              row.names = phenoData$sampleID)
mixedScore_scaled <- apply(mixedScore,2,function(x){x / max(x)})
svg("Plots/2_Cluster_heatmap.svg")
pheatmap(t(mixedScore_scaled),cluster_cols = fpcRes$FPC_res[[nClustInd]]$result$result,
         annotation_col = sampleAnnotation,show_rownames = FALSE,
         show_colnames = FALSE,clustering_method = "ward.D2",cutree_cols = nClusts[nClustInd])

dev.off(dev.cur())
png("Plots/2_Cluster_heatmap2.png",width = 35,height = 25,units = "cm",res = 300)
pheatmap(t(mixedScore_scaled),cluster_cols = fpcRes$FPC_res[[nClustInd]]$result$result,
         annotation_col = sampleAnnotation,show_rownames = FALSE,
         show_colnames = FALSE,clustering_method = "ward.D2",cutree_cols = nClusts[nClustInd])
dev.off(dev.cur())


# jpeg("Plots/2_Cluster_heatmap_tmp.jpeg",width = 1080,height = 1080,quality = 100)
# pheatmap(t(mixedScore),
#          annotation_col = sampleAnnotation,show_rownames = FALSE,
#          show_colnames = FALSE,clustering_method = "ward.D2",cutree_cols = nClusts[nClustInd],
#          clustering_distance_cols = "manhattan")
# dev.off(dev.cur())

##- (3) Cluster characterization -#
#-- Tissue purity -#

clusterTable <- data.frame(Clusters = as.character(unique(clusters)),
                           Tissue = rep("",length(unique(clusters))),
                           Purity = rep(0,length(unique(clusters))))
for (i in 1:nrow(clusterTable)){
  clust <- as.numeric(clusterTable$Clusters[i])
  clusterTable$Tissue[i] <- names(which.max(table(clusters == clust,phenoData$tissue)["TRUE",]))
  clusterTable$Purity[i] <-  round(sum(clusters == clust & phenoData$tissue == clusterTable$Tissue[i])/sum(clusters == clust),2)
}
clusterTable

#-- Tissue comparisons --#
plotData <- celScores |> 
  as.data.frame() |> 
  select(c(Task_1,Task_326,Task_275)) |> 
  mutate(Tissue = ifelse(phenoData$tissue == "Muscle","Skeletal Muscle",phenoData$tissue)) |> 
  mutate(Tissue = factor(Tissue,
                         levels = c("Skeletal Muscle", "Heart","Adipose Tissue",
                                    "Colon","Liver"))) 
margintop = 0
marginb = 3
g1 <- ggboxplot(data = plotData,x = "Tissue",y = "Task_1",fill = "grey",#fill = "Tissue",
                ylab = "Aerobic rephosphorylation of ATP \nfrom Glucose",size = 0.2) +
  theme(text = element_text(size = 10),
        plot.margin = margin(t = margintop,b = -marginb,unit = "cm")) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("Heart","Skeletal Muscle"),
                                        c("Heart","Adipose Tissue"),
                                        c("Heart","Colon"),
                                        c("Heart","Liver")),
                     size = 4,
                     vjust = 0.5,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  guides(x = guide_axis(angle = 45))

g1
ggsave("Plots/3_Heart_muscle_task.png",g1)

g2 <- ggboxplot(data = plotData,x = "Tissue",y = "Task_326",fill = "grey",#fill = "Tissue",
                ylab = "Glycerol-3-phosphate synthesis",size = 0.2) +
  theme(text = element_text(size = 10),
        plot.margin = margin(t = margintop,b = -marginb,unit = "cm")) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("Heart","Adipose Tissue"),
                                        c("Adipose Tissue","Skeletal Muscle"),
                                        c("Adipose Tissue","Colon"),
                                        c("Adipose Tissue","Liver")),
                     size = 4,
                     vjust = 0.5,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  theme(legend.position = "none") +
  guides(x = guide_axis(angle = 45))
g2
ggsave("Plots/3_Adipose_task.png",g2)

g3 <- ggboxplot(data = plotData,x = "Tissue",y = "Task_275",fill = "grey",#fill = "Tissue",
                ylab = "Fructose to glucose conversion",size = 0.2) +
  theme(text = element_text(size = 10),
        plot.margin = margin(t = margintop,b = -marginb,unit = "cm")) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("Liver","Colon"),
                                        c("Liver","Adipose Tissue"),
                                        c("Liver","Heart"),
                                        c("Liver","Skeletal Muscle")),
                     size = 4,
                     vjust = 0.5,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  theme(legend.position = "none") +
  guides(x = guide_axis(angle = 45))
g3
ggsave("Plots/3_Liver_task.png",g3)
g4 <- g1 + g2 + g3
ggsave("Plots/3_Tissue_tasks.png",g4,dpi = "retina",width = 16,height = 9,units = "cm")
# Heart and Muscle

# Adipose tissue

# Liver


#save("createdPlots","Plots/createdPlots.RData")
#save.image("../../Data/saved_wkspaces/1_end.RData")