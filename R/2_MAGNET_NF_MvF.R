#----------------------------------#
#   2- Comparing the task scores   #
#      between male/females in     #
#      the control MAGNET data     #
#----------------------------------#
#Last ran 11/09/2023
##- (1) Setup -#
set.seed(20230525)
#-- LocalDir definition
localDir <- "../Results/2_MAGNET_NF_MvF/"
setwd(localDir)
#-- Libraries -#
library(ggplot2)
library(tidyr)
library(purrr)
library("stringr")
library("pcaMethods")
library(svglite)
library(pcaMethods)
library("quantreg")
library(MASS)
library(glue)
library(ggprism)
source("../../R/R_subscript/task_analysis_load.R")


#-- Data -#
celScores <- t(read.table("Cellfie_scores/cellfie_score.txt",sep = ",",header = FALSE))
celScores_binary <- t(read.table("Cellfie_scores/cellfie_score_binary.txt",sep = ",",header = FALSE))
phenoData <- read.table("../0_Clean_Data/magnet_phenoData_NF.csv",sep = ",",header = T)
taskReport <- read.table("../0_Clean_Data/task_report.txt",
                         sep = ",",header = F)
taskReport <- taskReport[!taskReport$V1 %in% c("277","258"),]

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

mixedScore <- celScores * celScores_binary

onePct <- round(0.01*nrow(celScores_binary))
ubiquitous <- colSums(celScores_binary) <= onePct| colSums(celScores_binary) == nrow(celScores_binary) 


mixedScore_reduced <- mixedScore[,!ubiquitous]
celScores_reduced <- celScores[,!ubiquitous]
taskReport_reduced <- taskReport[!ubiquitous,]



##- (2) Data exploration -#
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


##- (3) hypothesis testing -#
test_gender <- wilcoxon_cliff(celScores,phenoData$gender,"Female",correction = "BH")
test_gender$Task_ID <- taskReport$Task_ID
test_gender$Task_name <- str_remove(taskReport$Task_name,"\\(.*$")
write.table(test_gender,"Tables/mann_whithney_M_v_F.txt",sep = "\t",col.names = T,row.names = F)

test_Hypertension <- wilcoxon_cliff(celScores,phenoData$Hypertension,"Yes",correction = "BH")
test_Hypertension$Task_ID <- taskReport$Task_ID
test_Hypertension$Task_name <- str_remove(taskReport$Task_name,"\\(.*$")
write.table(test_Hypertension,"Tables/mann_whithney_noHT_v_HT.txt",sep = "\t",col.names = T,row.names = F)

test_Race <- wilcoxon_cliff(celScores,phenoData$race,"AA",correction = "BH")
test_Race$Task_ID <- taskReport$Task_ID
test_Race$Task_name <- str_remove(taskReport$Task_name,"\\(.*$")
write.table(test_Race,"Tables/mann_whithney_Caucasian_v_AA.txt",sep = "\t",col.names = T,row.names = F)

test_Diabetes <- wilcoxon_cliff(celScores,phenoData$Diabetes,"Yes",correction = "BH")
test_Diabetes$Task_ID <- taskReport$Task_ID
test_Diabetes$Task_name <- str_remove(taskReport$Task_name,"\\(.*$")
write.table(test_Diabetes,"Tables/mann_whithney_noT2DM_v_T2DM.txt",sep = "\t",col.names = T,row.names = F)




lmData <- phenoData[,-c(1,5,6,7,9)]



test_regression <- quant_reg_clean(celScores,lmData)

write.table(test_regression,"Tables/quant_reg_all.txt",sep = "\t",col.names = T,row.names = F)


#-- Menopausal women
test_menopause <- wilcoxon_cliff(celScores[phenoData$gender == "Female",],phenoData[phenoData$gender == "Female",]$age > 45,TRUE)
test_menopause$Task_ID <- taskReport$Task_ID
test_menopause$Task_name <- str_remove(taskReport$Task_name,"\\(.*$")
write.table(test_menopause,"Tables/mann_whithney_under45_v_above45.txt",sep = "\t",col.names = T,row.names = F)

test_men_opause <- wilcoxon_cliff(celScores[phenoData$gender == "Male",],phenoData[phenoData$gender == "Male",]$age > 45,TRUE)
test_men_opause$Task_ID <- taskReport$Task_ID
test_men_opause$Task_name <- str_remove(taskReport$Task_name,"\\(.*$")
write.table(test_men_opause,"Tables/mann_whithney_under45_v_above45_males.txt",sep = "\t",col.names = T,row.names = F)


##- (4) Visualization

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


for(i in 1:nrow(test_menopause)){
  if(test_menopause$q.value[i] < 0.05){
    plotData <- data.frame(score = celScores[,i],
                           Sex = phenoData$gender,
                           Age = phenoData$age,
                           Age_45 = phenoData$age > 45)
    
    sigTable <- data.frame(
      group1 = c(FALSE,FALSE),
      group2 = c(TRUE,TRUE),
      p.adj = c(test_menopause$q.value[i],test_men_opause$q.value[i]),
      y.position = c(max(celScores[,i])*1.05,max(celScores[,i])*1.05),
      Sex =  c("Female","Male"),
      Age_45 = c(FALSE,TRUE)
    )
    ggplot(plotData,aes(Age,score,color = Sex))+
      geom_point()+
      geom_vline(xintercept = 45,alpha = 0.5)+
      facet_wrap(vars(Sex))+
      ylab("CellFie Score")+
      ggtitle(paste0(colnames(celScores)[i],": ",formatedTaskName[i]))+
      theme_cardio()
    ggsave(paste0("Plots/menopause/menopause_investigation_point_",colnames(celScores)[i],".svg"))
    
    ggplot(plotData,aes(Age_45,score,fill = Age_45))+
      geom_boxplot()+
      geom_violin(linewidth = 0.2,fill = NA)+
      #geom_point(color = "black",alpha = 0.4,position = position_jitter(height = 0,width = 0.05))+

      facet_wrap(vars(Sex))+
      ylab("CellFie Score")+
      ggtitle(paste0(colnames(celScores)[i],": ",formatedTaskName[i]))+
      theme_cardio()+
      labs(fill = "Older than 45")+
      xlab("Older than 45")+
      add_pvalue(sigTable,label = "FDR = {format(p.adj,digits = 3,scientific = -2)}",tip.length = 0.01)
    ggsave(paste0("Plots/menopause/menopause_investigation_box_",colnames(celScores)[i],".svg"))
    }
}

#-- Race

for(i in 1:nrow(test_Race)){
  if(test_Race$q.value[i] < 0.05){
    
    plotData <- data.frame(score = celScores[,i],
                           Sex = phenoData$gender,
                           Age = phenoData$age,
                           Race = phenoData$race,
                           Age_45 = phenoData$age > 45 )
    
    sigTable <- data.frame(
      group1 = "Caucasian",
      group2 = "AA",
      p.adj = test_Race$q.value[i],
      y.position = c(max(celScores[,i])*1.05)
    )
    ggplot(plotData,aes(Race,score))+
      geom_boxplot(aes(fill = Race))+
      #geom_point(color = "black",alpha = 0.4,position = position_jitter(height = 0,width = 0.05))+
      geom_violin(linewidth = 0.2,fill = NA)+
      ylab("CellFie Score")+
      ggtitle(paste0(colnames(celScores)[i],": ",formatedTaskName[i]))+
      theme_cardio()+
      labs(fill = "Ethnicity")+
      xlab("Ethnicity")+
      add_pvalue(data = sigTable,"FDR = {format(p.adj,digits = 3,scientific = -2)}",tip.length = 0.01)
    ggsave(paste0("Plots/ethnicity/ethnicity_investigation_box_",colnames(celScores)[i],".svg"))
  }
}

#-- gender

for(i in 1:nrow(test_gender)){
  if(test_gender$q.value[i] < 0.05){
    plotData <- data.frame(score = celScores[,i],
                           Sex = phenoData$gender,
                           Age = phenoData$age,
                           Race = phenoData$race,
                           Age_45 = phenoData$age > 45 )
    sigTable <- data.frame(
      group1 = "Male",
      group2 = "Female",
      p.adj = test_gender$q.value[i],
      y.position = c(max(celScores[,i])*1.05)
    )
    
    ggplot(plotData,aes(Sex,score))+
      geom_boxplot(aes(fill = Sex))+
      geom_violin(linewidth = 0.2,fill = NA)+
      ylab("CellFie Score")+
      ggtitle(paste0(colnames(celScores)[i],": ",formatedTaskName[i]))+
      theme_cardio()+
      labs(fill = "Sex")+
      xlab("Sex")+
      add_pvalue(data = sigTable,"FDR = {format(p.adj,digits = 3,scientific = -2)}",tip.length = 0.01)
    ggsave(paste0("Plots/sex/sex_investigation_box_",colnames(celScores)[i],".svg"))
  }
}

save.image("../../Data/saved_wkspaces/2_end")
