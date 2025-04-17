#----------------------------------#
#   2- Comparing the task scores   #
#      between DCM/Control in      #
#      the  MAGNET data            #
#----------------------------------#

##- (1) Setup -#
set.seed(20230525)
#-- LocalDir definition
localDir <- "../Results/2_MAGNET_NFvDCM/"
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
library("pcaMethods")
library(svglite)
library("quantreg")
library(MASS)
library(glue)
library(ggprism)
library(fpc)
library(gtsummary)
library(ggpubr)
source("../../R/R_subscript/task_analysis_load.R")


#-- Data -#
celScores <- t(read.table("Cellfie_scores/cellfie_score.txt",sep = ",",header = FALSE))
celScores_binary <- t(read.table("Cellfie_scores/cellfie_score_binary.txt",sep = ",",header = FALSE))
phenoData <- read.table("../0_Clean_Data/magnet_phenoData_complete.csv",sep = ",",header = T)
taskReport <- read.table("../0_Clean_Data/task_report.txt",
                         sep = ",",header = F)
taskReport <- taskReport[!taskReport$V1 %in% c("277","258"),]

TOI <- read.xlsx("../../Data/MetabolicTasks_closer look choice.xlsx",sheet = 1)
TOI <- TOI[!TOI$ID %in% c("277","258"),]
TOI[TOI$ID == "252",]$Task.of.interest <- TRUE 


#-- Data reformatting -#
# when no information on HT or diabetes, assumed non-affected
phenoData$Hypertension[is.na(phenoData$Hypertension)] <- "No"
phenoData$Diabetes[is.na(phenoData$Diabetes)] <- "No"

# reformatting task repport
colnames(taskReport) <- c("Task_ID","Subsystem","Source","Task_name","Failed")
taskReport <- taskReport[1:(nrow(taskReport)-2),]
taskReport$Task_ID <- paste0("Task_",taskReport$Task_ID)
taskReport$Failed <- sapply(taskReport$Failed,function(x){
  switch(x,
         "true" = F,
         T)
})

# Reformatting scores
colnames(celScores) <- colnames(celScores_binary) <- taskReport$Task_ID
rownames(celScores) <- rownames(celScores_binary) <- phenoData$sampleID

# Failed tasks
minusOneScore <- celScores_binary[1,] == -1
cat(paste0(sum(taskReport$Failed)," tasks do not have score due to having failed, \n",
           sum(minusOneScore) - sum(taskReport$Failed)," tasks do not have score for other reasons, ",
           "such as not having essential rxns with associated GPR"))

# Removing failed tasks
celScores <- celScores[,!minusOneScore]
celScores_binary <- celScores_binary[,!minusOneScore]
TOI <- TOI[!minusOneScore,]
taskReport <- taskReport[!minusOneScore,]

# Generating mixed score
mixedScore <- celScores * celScores_binary

# Defining ubiquitous tasks
onePct <- round(0.01*nrow(celScores_binary))
ubiquitous <- colSums(celScores_binary) <= onePct| colSums(celScores_binary) == nrow(celScores_binary)

# Scores with ubiquitous tasks removed
mixedScore_reduced <- mixedScore[,!ubiquitous]
celScores_reduced <- celScores[,!ubiquitous]
taskReport_reduced <- taskReport[!ubiquitous,]



##- (2) Data exploration -#
#-- Subject information
summary_table <- phenoData |> 
  gtsummary::tbl_summary(by = etiology,
                         include = c("age","gender","race",
                                     "Diabetes","Hypertension")) |> 
  add_p()
summary_table
gt::gtsave(as_gt(summary_table),
           "Tables/Summary_subjects.docx")
#-- Preparing the PCA plot Data


pcaRes <- pcaMethods::pca(mixedScore_reduced)

plotData <- data.frame(PC1 = pcaRes@scores[,1],
                       PC2 = pcaRes@scores[,2])
plotData <- cbind(plotData,phenoData)

#-- Plotting
#Gender
ggplot(plotData, aes(PC1,PC2,color = gender))+
  geom_point()+
  xlab(paste0("PC1 (",round(pcaRes@R2[1],1)*100,"%)"))+
  ylab(paste0("PC2 (",round(pcaRes@R2[2],1)*100,"%)"))+
  theme_cardio()
#Ethnicity
ggplot(plotData, aes(PC1,PC2,color = race))+
  geom_point()+
  xlab(paste0("PC1 (",round(pcaRes@R2[1],1)*100,"%)"))+
  ylab(paste0("PC2 (",round(pcaRes@R2[2],1)*100,"%)"))+
  theme_cardio()
#Age
ggplot(plotData, aes(PC1,PC2,color = age))+
  geom_point()+
  xlab(paste0("PC1 (",round(pcaRes@R2[1],1)*100,"%)"))+
  ylab(paste0("PC2 (",round(pcaRes@R2[2],1)*100,"%)"))+
  theme_cardio()
#Hypertension
ggplot(plotData, aes(PC1,PC2,color = Hypertension))+
  geom_point()+
  xlab(paste0("PC1 (",round(pcaRes@R2[1],1)*100,"%)"))+
  ylab(paste0("PC2 (",round(pcaRes@R2[2],1)*100,"%)"))+
  theme_cardio()

#Diabetes
ggplot(plotData, aes(PC1,PC2,color = Diabetes))+
  geom_point()+
  xlab(paste0("PC1 (",round(pcaRes@R2[1],1)*100,"%)"))+
  ylab(paste0("PC2 (",round(pcaRes@R2[2],1)*100,"%)"))+
  theme_cardio()

#Weight
ggplot(plotData, aes(PC1,PC2,color = weight))+
  geom_point()+
  xlab(paste0("PC1 (",round(pcaRes@R2[1],1)*100,"%)"))+
  ylab(paste0("PC2 (",round(pcaRes@R2[2],1)*100,"%)"))+
  theme_cardio()
#BMI
ggplot(plotData, aes(PC1,PC2,color = BMI))+
  geom_point()+
  xlab(paste0("PC1 (",round(pcaRes@R2[1],1)*100,"%)"))+
  ylab(paste0("PC2 (",round(pcaRes@R2[2],1)*100,"%)"))+
  theme_cardio()
#Disease
ggplot(plotData, aes(PC1,PC2,color = etiology))+
  geom_point()+
  xlab(paste0("PC1 (",round(pcaRes@R2[1],1)*100,"%)"))+
  ylab(paste0("PC2 (",round(pcaRes@R2[2],1)*100,"%)"))+
  theme_cardio()

# -> we observe that even in mixed score space, the two first PCs are dominated by patient's etiology


##- (3) Hypothesis testing -#
#Fixing mixed scores to avoid ties
# binary2 <- celScores_binary[,!ubiquitous]
# binary2[binary2 == 0] <- 0.0001
# mixedScore2 <- celScores * binary2

test_etiology <- wilcoxon_cliff(celScores,phenoData$etiology,"DCM",correction = "BH")
test_etiology$Task_ID <- taskReport$Task_ID
test_etiology$Task_name <- str_remove(taskReport$Task_name,"\\(.*$")
test_etiology$TOI <- TOI$Task.of.interest
test_etiology$In_Clustering <- !ubiquitous
# test_etiology <- wilcoxon_cliff(celScores,phenoData$etiology,"DCM",correction = "BH")
# test_etiology$Task_ID <- taskReport$Task_ID
# test_etiology$Task_name <- str_remove(taskReport$Task_name,"\\(.*$")
write.table(test_etiology,"Tables/mann_whithney_NF_v_DCM.txt",sep = "\t",col.names = T,row.names = F)
# 
# #What if linear
# test_etiology_lin <- wilcoxon_cliff(celScores,phenoData$etiology,"DCM",correction = "BH")
# test_etiology_lin$Task_ID <- taskReport$Task_ID
# test_etiology_lin$Task_name <- str_remove(taskReport$Task_name,"\\(.*$")
# 
# table(test_etiology$q.value < 0.05, test_etiology_lin$q.value < 0.05)
##- (4 tentative) Clustering DCM scores -#
#-- PCA -#
DCMScores <- mixedScore[phenoData$etiology == "DCM",]
pcaRes <- pcaMethods::pca(DCMScores,nPcs = ncol(DCMScores))
pcaScores <- pcaRes@scores[,pcaRes@R2 > 0.01]

#-- FPC run -#
B <- 100000 #Number of boostraps
nClusts <- 2:15 #Number of clusters to attempt

fpcRes <- fpc_diagnostic_parallel(dist(pcaScores),B,nClusts,seed = 20222023)

#-- FPC Visualization -#
p1 <- ggplot(fpcRes$plotData,aes(x=nClust,y=Mean_jaccard))+
  geom_point()+
  geom_hline(yintercept = 0.6)+
  theme_cardio()+
  ylab("Mean JC across bootstraps")+
  xlab("Number of clusters")

plot(p1)
ggsave("Plots/1_FPC_plot.svg",p1)

#-- Clustering -#
allAbove <- rep(F,length(nClusts))
for (i in 1:length(nClusts)){
  n <- nClusts[i]
  allAbove[i] <- sum(fpcRes$plotData$nClust == n & fpcRes$plotData$Mean_jaccard > 0.6) == n
}
if(sum(allAbove) == 0){
  cat("There are no stable clusters")
}else{
  nClustInd <- max(which(allAbove))
  
  cat(paste0("There are ",nClusts[nClustInd]," stable clusters"))
  
  clusters <- fpcRes$FPC_res[[nClustInd]]$partition
  clustTable <- data.frame(Subject = phenoData[phenoData$etiology == "DCM",]$ID,
                           Cluster = clusters)
  write.table(clustTable,"Tables/DCM_clusts_from_whole.txt",sep = "\t",col.names = T,row.names = F)
}



##- (4) Vizualization

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
formatedTaskIDs <- taskReport$Task_ID
formatedTaskIDs <- str_remove(formatedTaskIDs,"Task_")
formatedTaskIDs <- str_pad(formatedTaskIDs,3,"left","0")
formatedTaskIDs <- paste0("Task_",formatedTaskIDs)
for(i in 1:nrow(test_etiology)){
  if(test_etiology$q.value[i] < 0.05){
    plotData <- data.frame(score = celScores[,i],
                           Sex = phenoData$gender,
                           Age = phenoData$age,
                           Race = phenoData$race,
                           Group = phenoData$etiology)
    plotData$Group[plotData$Group == "NF"] <- "Control"
    sigTable <- data.frame(
      group1 = "Control",
      group2 = "DCM",
      p.adj = test_etiology$q.value[i],
      y.position = c(max(celScores[,i])*1.05)
    )
    ggplot(plotData,aes(Group,score))+
      geom_boxplot(aes(fill = Group))+
      geom_violin(linewidth = 0.2,fill = NA)+
      #geom_point(color = "black",alpha = 0.4,position = position_jitter(height = 0,width = 0.05))+
      ylab("CellFie Score")+
      ggtitle(paste0(colnames(celScores)[i],": ",formatedTaskName[i]))+
      theme_cardio()+
      labs(fill = "Group")+
      xlab("Group")+
      add_pvalue(data = sigTable,"FDR = {format(p.adj,digits = 3,scientific = -2)}",tip.length = 0.01)
    
    ggsave(paste0("Plots/mann_whitney/",formatedTaskIDs[i],"_donor_v_DCM.svg"),height = 10,width = 12,units = "cm")
  }
}

ggplot(test_etiology,aes(logFC,-log(p.value))) +
  geom_point() +
  geom_hline(yintercept = -log(0.05),color = "red")
ggsave("Plots/volcano_plot_pvals.png")

ggplot(test_etiology,aes(logFC,-log(q.value))) +
  geom_point() +
  geom_hline(yintercept = -log(0.05),color = "red")
ggsave("Plots/volcano_plot_qvals.png")

######################"
# test_etiology$scaled.median.diff <-test_etiology$median.diff / apply(mixedScore2,2,max)
# g1 <- ggplot(test_etiology,aes(cliff.delta,-log10(p.value)))+
#   geom_point()+
#   geom_hline(yintercept = -log10(0.05),color = "red")
# g1
# g2 <- ggplot(test_etiology,aes(median.diff,-log10(p.value)))+
#   geom_point()+
#   geom_hline(yintercept = -log10(0.05),color = "red")
# g2
# g3 <- ggplot(test_etiology,aes(cliff.delta,median.diff))+
#   geom_point()
# g3
# g4 <- ggplot(test_etiology,aes(cliff.delta,scaled.median.diff))+
#   geom_point()
# g4
# g5 <- ggplot(test_etiology,aes(scaled.median.diff,-log10(p.value)))+
#   geom_point()+
#   geom_hline(yintercept = -log10(0.05),color = "red")
# g5
# 
# library(patchwork)
# (g1 + g2)/g3
# g4 <- (g1 + g2)/g3
# ggsave("Plots/sigplot.svg",g4)


save.image("../../Data/saved_wkspaces/2_end.RData")
#load("../../Data/saved_wkspaces/2_end.RData")

#------------------------

celScores[phenoData$etiology == "NF",] |> 
  as.data.frame() |> 
  map_dbl(\(x) {
    x <- shapiro.test(x)
    return(x$p.value)
    }) |>
  data.frame(`p.value` = _) |> 
  ggplot(aes(-log(`p.value`))) +
  geom_histogram() +
  geom_vline(xintercept = -log(0.05)) +
  theme_pubr()


dat <- t(mixedScore_reduced) 
colnames(dat) <- phenoData$ID
annotData <- data.frame(row.names = colnames(dat),
                        Etiology = phenoData$etiology)
pheatmap::pheatmap(dat,clustering_method = "ward.D2",
                   annotation_col = annotData)
  
