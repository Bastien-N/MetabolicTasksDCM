#------------------------------------#
#   6- Assembling results together   #
#------------------------------------#
##- (1) Setup -#
set.seed(20230525)
#-- LocalDir definition
localDir <- "../Results/6_Final_touches/"
setwd(localDir)
#-- Libraries -#
library(ggplot2)
library(tidyr)
library(purrr)
library("stringr")
library(svglite)
library(pcaMethods)
library(glue)
library(ggprism)
library(openxlsx)
library(dplyr)
library(ggpubr)
library(ggbreak)
library(patchwork)
library(ggtext)

#-- Data -#
prelim_lit <- openxlsx::read.xlsx("../../Data/Prelim_task_literature_link.xlsx")
taskReport <- read.table("../0_Clean_Data/task_report.txt",
                         sep = ",",header = F)

taskReport <- taskReport[!taskReport$V1 %in% c("277","258"),]


magnet_NF_sex <- read.table("../2_MAGNEt_NF_MvF/Tables/mann_whithney_M_v_F.txt",sep = "\t",header = T)
magnet_NF_ethnicity <- read.table("../2_MAGNEt_NF_MvF/Tables/mann_whithney_Caucasian_v_AA.txt",sep = "\t",header = T) 

magnet_DCM <- read.table("../3_MAGNET_NFvDCM/Tables/mann_whithney_NF_v_DCM.txt",sep = "\t",header = T)
magnet_C1_v_C2 <- read.table("../4_MAGNET_DCM_clustering/Tables/mann_whitney_cluster_tasks.txt",sep = "\t",header = T)
magnet_Control_v_Clust <- read.table("../5_MAGNET_clusts_v_NF/Tables/mann_whitney_tasks_clusters_v_control.txt",sep = "\t",header = T)

magnet_NF_celScore_binary <- t(read.table("../2_MAGNEt_NF_MvF/Cellfie_scores/cellfie_score_binary.txt",sep = ",",header = FALSE))
magnet_all_celScore_binary <- t(read.table("../3_MAGNEt_NFvDCM/Cellfie_scores/cellfie_score_binary.txt",sep = ",",header = FALSE))
magnet_DCM_celScore_binary <- t(read.table("../4_MAGNEt_DCM_clustering/Cellfie_scores/cellfie_score_binary.txt",sep = ",",header = FALSE))

TOI <- read.xlsx("../../Data/MetabolicTasks_closer look choice.xlsx",sheet = 1)
TOI <- TOI[!TOI$ID %in% c("277","258"),]
TOI[TOI$ID == "252",]$Task.of.interest <- TRUE 
#-- Data reformatting -#
taskReport <- taskReport[-((nrow(taskReport)-1):nrow(taskReport)),]
colnames(taskReport) <- c("Task_ID","Subsystem","Source","Task_name","Should_fail")
taskReport$Task_ID <- paste0("Task_",taskReport$Task_ID)

taskReport$Subsystem <- str_replace(taskReport$Subsystem,"\\&","and")

existingSubs <- unique(taskReport$Subsystem)

taskReport$Task_name <- str_replace(taskReport$Task_name,
                                    "\\(complete oxidation\\)",
                                    "complete oxidation")
taskReport$Task_name <- str_remove(taskReport$Task_name,"\\(.*$")




magnet_NF_pctActive <- colSums(magnet_NF_celScore_binary == 1) / nrow(magnet_NF_celScore_binary)
magnet_all_pctActive <- colSums(magnet_all_celScore_binary == 1) / nrow(magnet_all_celScore_binary)
magnet_DCM_pctActive <- colSums(magnet_DCM_celScore_binary == 1) / nrow(magnet_DCM_celScore_binary)

##- (2) Contructing tables

#-- Result table 1 (Magnet NF)

ResultsTable1 <- data.frame(Task_ID = taskReport$Task_ID,
                            Task_Name = taskReport$Task_name,
                            Subsystem = taskReport$Subsystem,
                            Additional_note = NA,
                            Potential_DCM_link = NA,
                            Pct_active_Controls = magnet_NF_pctActive,
                            M_v_F_FDR = NA,
                            M_v_F_Cliff = NA,
                            M_v_F_FC = NA,
                            C_v_AA_FDR = NA,
                            C_v_AA_Cliff = NA,
                            C_v_AA_FC = NA)
format1 <- c(NA,NA,"byType",NA,NA,NA,"fdr","Cliff","FC","fdr","Cliff","FC")
#Prelim literature search
for (i in 1:nrow(prelim_lit)){
  if(!is.na(prelim_lit$Added.context[i])){
    ResultsTable1$Additional_note[i] <- prelim_lit$Added.context[i]
  }
  if(!is.na(prelim_lit$Possible.DCM.connection[i])){
    ResultsTable1$Potential_DCM_link[i] <- prelim_lit$Possible.DCM.connection[i]
  }
}

#Numbers
for(i in 1:nrow(magnet_NF_sex)){
  ind <- which(ResultsTable1$Task_ID == magnet_NF_sex$Task_ID[i])
  ResultsTable1$M_v_F_FDR[ind] <- magnet_NF_sex$q.value[i]
  ResultsTable1$M_v_F_Cliff[ind] <- magnet_NF_sex$cliff.delta[i]
  ResultsTable1$M_v_F_FC[ind] <- magnet_NF_sex$logFC[i] |> 
    (function(x){2^x})() |> 
    (function(x){ifelse(x < 1,-1/x,x)})() 
  
  ResultsTable1$C_v_AA_FDR[ind] <- magnet_NF_ethnicity$q.value[i]
  ResultsTable1$C_v_AA_Cliff[ind] <- magnet_NF_ethnicity$cliff.delta[i]
  ResultsTable1$C_v_AA_FC[ind] <- magnet_NF_ethnicity$logFC[i] |> 
    (function(x){2^x})() |> 
    (function(x){ifelse(x < 1,-1/x,x)})()
}

#-- Results table 2 (Magnet NF_v_DCM and clusts)

ResultsTable2 <- data.frame(Task_ID = taskReport$Task_ID,
                            Task_Name = taskReport$Task_name,
                            Subsystem = taskReport$Subsystem,
                            Additional_note = NA,
                            Potential_DCM_link = NA,
                            Pct_active_all = magnet_all_pctActive,
                            Control_v_DCM_FDR = NA,
                            Control_v_DCM_Cliff = NA,
                            Control_v_DCM_FC = NA,
                            Control_v_C1_FDR = NA,
                            Control_v_C1_Cliff = NA,
                            Control_v_C1_FC = NA,
                            Control_v_C2_FDR = NA,
                            Control_v_C2_Cliff = NA,
                            Control_v_C2_FC = NA,
                            Pct_active_DCM = magnet_DCM_pctActive,
                            C1_v_C2_FDR = NA,
                            C1_v_C2_Cliff = NA,
                            C1_v_C2_FC = NA,
                            Sig_description = NA,
                            TOI = TOI$Task.of.interest)

format2 <- c(NA,NA,"byType",NA,NA,
             NA,
             "fdr","Cliff","FC",
             "fdr","Cliff","FC",
             "fdr","Cliff","FC",NA,
             "fdr","Cliff","FC",NA,NA)

#Numbers
for(i in 1:nrow(magnet_DCM)){
  ind <- which(ResultsTable2$Task_ID == magnet_DCM$Task_ID[i])
  ResultsTable2$Control_v_DCM_FDR[ind] <- magnet_DCM$q.value[i]
  ResultsTable2$Control_v_DCM_Cliff[ind] <- magnet_DCM$cliff.delta[i]
  ResultsTable2$Control_v_DCM_FC[ind] <- magnet_DCM$logFC[i] |> 
    (function(x){2^x})() |> 
    (function(x){ifelse(x < 1,-1/x,x)})()
}
for(i in 1:nrow(magnet_C1_v_C2)){
  ind <- which(ResultsTable2$Task_ID == magnet_C1_v_C2$Task_ID[i])
  ResultsTable2$C1_v_C2_FDR[ind] <- magnet_C1_v_C2$q.value[i]
  ResultsTable2$C1_v_C2_Cliff[ind] <- magnet_C1_v_C2$cliff.delta[i]
  ResultsTable2$C1_v_C2_FC[ind] <- magnet_C1_v_C2$logFC[i] |> 
    (function(x){2^x})() |> 
    (function(x){ifelse(x < 1,-1/x,x)})()
}
for(i in 1:nrow(magnet_Control_v_Clust)){
  ind <- which(ResultsTable2$Task_ID == magnet_Control_v_Clust$Task_ID[i])
  if(magnet_Control_v_Clust$Comp[i] == "Control_vs_Cluster1"){
    ResultsTable2$Control_v_C1_FDR[ind] <- magnet_Control_v_Clust$q.value[i]
    ResultsTable2$Control_v_C1_Cliff[ind] <- magnet_Control_v_Clust$cliff.delta[i]
    ResultsTable2$Control_v_C1_FC[ind] <- magnet_Control_v_Clust$logFC[i] |> 
      (function(x){2^x})() |> 
      (function(x){ifelse(x < 1,-1/x,x)})()
  }else{
    ResultsTable2$Control_v_C2_FDR[ind] <- magnet_Control_v_Clust$q.value[i]
    ResultsTable2$Control_v_C2_Cliff[ind] <- magnet_Control_v_Clust$cliff.delta[i]
    ResultsTable2$Control_v_C2_FC[ind] <- magnet_Control_v_Clust$logFC[i] |> 
      (function(x){2^x})() |> 
      (function(x){ifelse(x < 1,-1/x,x)})()
  }

}

#-- Examining the relationships
sig_overall <- !is.na(ResultsTable2$Control_v_DCM_FDR) & ResultsTable2$Control_v_DCM_FDR < 0.05
sig_C1 <- !is.na(ResultsTable2$Control_v_C1_FDR) & ResultsTable2$Control_v_C1_FDR < 0.05
sig_C2 <- !is.na(ResultsTable2$Control_v_C2_FDR) & ResultsTable2$Control_v_C2_FDR < 0.05
sig_between <- !is.na(ResultsTable2$C1_v_C2_FDR) & ResultsTable2$C1_v_C2_FDR < 0.05
oppositeSigns <- !is.na(ResultsTable2$Control_v_C1_Cliff) & sign(ResultsTable2$Control_v_C1_Cliff) != sign(ResultsTable2$Control_v_C2_Cliff)

inAgreement <- sig_overall & sig_C1 & sig_C2 & !oppositeSigns
oppositeSig <- sig_C1 & sig_C2 & oppositeSigns

unequalSig <- sig_overall & sig_between
C1sigDom <- sig_overall & sig_C1 & !sig_C2
C2sigDom <- sig_overall & sig_C2 & !sig_C1

onlyC1 <- !sig_overall & sig_C1 & !sig_C2
onlyC2 <- !sig_overall & sig_C2 & !sig_C1
sigBeforeClust <- sig_overall & !sig_C1 & !sig_C2

sigDescription <- rep(NA,nrow(ResultsTable2))
for(i in 1:length(sigDescription)){
  if(inAgreement[i]){
    sigDescription[i] <- "Significant in both Clusters"
  }else if(C1sigDom[i]){
    sigDescription[i] <- "Significant only because of Cluster 1"
  }else if(C2sigDom[i]){
    sigDescription[i] <- "Significant only because of Cluster 2"
  }else if(onlyC1[i]){
    sigDescription[i] <- "Only significant in Cluster 1"
  }else if(onlyC2[i]){
    sigDescription[i] <- "Only significant in Cluster 2"
  }else if(oppositeSig[i]){
    sigDescription[i] <- "Significant in opposite direction between clusters"
  }else if(sigBeforeClust[i]){
    sigDescription[i] <- "Only significant without clustering"
  }
}

ResultsTable2$Sig_description <- sigDescription


##- (3) Assembling excel file
#--Styles
baseStyle <- createStyle(valign = "top",wrapText = T)
SigQvalStyle <- createStyle(textDecoration = "bold",fontColour = "#ff0000")
SigCliffStyle <- createStyle(textDecoration = "bold")


subSystemColors <- c("#ccccff","#fa8072","#c1e1c3","#ffffcc","#dddddd","#facf5e",
                              "#98c1d9","#8ef0e7","#ec8c4f","#f683b2")

#-- wb first sheet
wb <- createWorkbook()
addWorksheet(wb,sheetName = "Magnet_step2")
writeDataTable(wb,sheet = "Magnet_step2",ResultsTable1,tableStyle = "none")
colWidth <- c(8.42,
              8.43*3,
              8.83*1.5,
              8.43*5,
              8.43*5)
for(i in 1:length(format1)){
  addStyle(wb,sheet = "Magnet_step2",style = baseStyle,rows = 1:353,cols = 1:10,gridExpand = T)
  switch (format1[i],
    "fdr" = conditionalFormatting(wb,"Magnet_step2",i,2:353,"< 0.05",style = SigQvalStyle,type = "expression"),
    "byType" = {
      for(ii in 1:length(subSystemColors)){
        tmpStyle <- createStyle(bgFill = subSystemColors[ii])
        conditionalFormatting(wb,"Magnet_step2",i,2:353,rule = existingSubs[ii],
                              style = tmpStyle,
                              type = "contains")
      }
      },
    "Cliff" = {
      tmpRule <- paste0("$",LETTERS[i-1],"2:$",LETTERS[i-1],"353 < 0.05")
      conditionalFormatting(wb,"Magnet_step2",i,2:353,tmpRule,style = SigCliffStyle,type = "expression")
      conditionalFormatting(wb,"Magnet_step2",i,2:353,rule = c(-1,0,1),
                            style = c("#0000ff","#ffffff","#cc0000"),type = "colourScale")
      },
    "FC" = {
      tmpRule <- paste0("$",LETTERS[i-2],"2:$",LETTERS[i-2],"353 < 0.05")
      extreme <- max(abs(ResultsTable1[,i]),na.rm = T)
      conditionalFormatting(wb,"Magnet_step2",i,2:353,tmpRule,style = SigCliffStyle,type = "expression")
      conditionalFormatting(wb,"Magnet_step2",i,2:353,rule = c(-extreme,0,extreme),
                            style = c("#0000ff","#ffffff","#cc0000"),type = "colourScale")
    }
    )
}
setColWidths(wb,sheet = "Magnet_step2",cols = 1:5,widths = colWidth)
setColWidths(wb,sheet = "Magnet_step2",cols = 6:ncol(ResultsTable1),widths = "auto")
addStyle(wb,sheet = "Magnet_step2",style = baseStyle,rows = 1:353,cols = 1:10,gridExpand = T)
freezePane(wb,sheet = "Magnet_step2",firstActiveCol = 2)
#-- wb second sheet 
addWorksheet(wb,sheetName = "Magnet_step3to5")
writeDataTable(wb,sheet = "Magnet_step3to5",ResultsTable2,tableStyle = "none")
colWidth <- c(8.42,
              8.43*3,
              8.83*1.5,
              8.43*5,
              8.43*5)
for(i in 1:length(format2)) {
  addStyle(wb,sheet = "Magnet_step3to5",style = baseStyle,rows = 1:353,cols = 1:10,gridExpand = T)
  switch (format2[i],
          "fdr" = conditionalFormatting(wb,"Magnet_step3to5",i,2:353,"< 0.05",style = SigQvalStyle,type = "expression"),
          "byType" = {
            for(ii in 1:length(subSystemColors)){
              tmpStyle <- createStyle(bgFill = subSystemColors[ii])
              conditionalFormatting(wb,"Magnet_step3to5",i,2:353,rule = existingSubs[ii],
                                    style = tmpStyle,
                                    type = "contains")
            }
          },
          "Cliff" = {
            tmpRule <- paste0("$",LETTERS[i-1],"2:$",LETTERS[i-1],"353 < 0.05")
            conditionalFormatting(wb,"Magnet_step3to5",i,2:353,tmpRule,style = SigCliffStyle,type = "expression")
            conditionalFormatting(wb,"Magnet_step3to5",i,2:353,rule = c(-1,0,1),
                                  style = c("#0000ff","#ffffff","#cc0000"),type = "colourScale")
          },
          "med.diff" = {
            tmpRule <- paste0("$",LETTERS[i-2],"2:$",LETTERS[i-2],"353 < 0.05")
            extreme <- max(abs(ResultsTable2[,i]),na.rm = T)
            conditionalFormatting(wb,"Magnet_step3to5",i,2:353,tmpRule,style = SigCliffStyle,type = "expression")
            conditionalFormatting(wb,"Magnet_step2",i,2:353,rule = c(-extreme,0,extreme),
                                  style = c("#0000ff","#ffffff","#cc0000"),type = "colourScale")
          }
  )
}
setColWidths(wb,sheet = "Magnet_step3to5",cols = 1:5,widths = colWidth)
setColWidths(wb,sheet = "Magnet_step3to5",cols = 6:ncol(ResultsTable2),widths = "auto")
addStyle(wb,sheet = "Magnet_step3to5",style = baseStyle,rows = 1:353,cols = 1:ncol(ResultsTable2),gridExpand = T)
freezePane(wb,sheet = "Magnet_step3to5",firstActiveCol = 2) 

saveWorkbook(wb,"agreggated_results_Tasks_v_groups.xlsx",overwrite = T)

# Second results table
tb <- data.frame(`Task` = paste0(ResultsTable2$Task_ID,": ",ResultsTable2$Task_Name),
                 `FC DCM vs Controls` = ResultsTable2$Control_v_DCM_FC,
                 `FC Cluster 1 vs Controls` = ResultsTable2$Control_v_C1_FC,
                 `FC Cluster 2 vs Controls` = ResultsTable2$Control_v_C2_FC,
                 `FC Cluster 2 vs Cluster 1` = ResultsTable2$C1_v_C2_FC)
tobold <- tb[,-1] > 1.1 | tb[,-1] < -1.1
tb$FC.DCM.vs.Controls <- format(tb$FC.DCM.vs.Controls,digits = 2,nsmall = 2)
tb$FC.Cluster.1.vs.Controls <- format(tb$FC.Cluster.1.vs.Controls,digits = 2,nsmall = 2)
tb$FC.Cluster.2.vs.Controls <- format(tb$FC.Cluster.2.vs.Controls,digits = 2,nsmall = 2)
tb$FC.Cluster.2.vs.Cluster.1 <- format(tb$FC.Cluster.2.vs.Cluster.1,digits = 2,nsmall = 2)




for ( i in 1:nrow(tb)) {
  tb[i,2] <- case_when(
    ResultsTable2$Control_v_DCM_FDR[i] < 0.001 ~ paste0(tb[i,2],"***"),
    ResultsTable2$Control_v_DCM_FDR[i] < 0.01 ~ paste0(tb[i,2],"**"),
    ResultsTable2$Control_v_DCM_FDR[i] < 0.05 ~ paste0(tb[i,2],"*"),
    .default = paste0(tb[i,2])
  )
  tb[i,3] <- case_when(
    ResultsTable2$Control_v_C1_FDR[i] < 0.001 ~ paste0(tb[i,3],"***"),
    ResultsTable2$Control_v_C1_FDR[i] < 0.01 ~ paste0(tb[i,3],"**"),
    ResultsTable2$Control_v_C1_FDR[i] < 0.05 ~ paste0(tb[i,3],"*"),
    .default = paste0(tb[i,3])
  )
  tb[i,4] <- case_when(
    ResultsTable2$Control_v_C2_FDR[i] < 0.001 ~ paste0(tb[i,4],"***"),
    ResultsTable2$Control_v_C2_FDR[i] < 0.01 ~ paste0(tb[i,4],"**"),
    ResultsTable2$Control_v_C2_FDR[i] < 0.05 ~ paste0(tb[i,4],"*"),
    .default = paste0(tb[i,4])
  )
  tb[i,5] <- case_when(
    ResultsTable2$C1_v_C2_FDR[i] < 0.001 ~ paste0(tb[i,5],"***"),
    ResultsTable2$C1_v_C2_FDR[i] < 0.01 ~ paste0(tb[i,5],"**"),
    ResultsTable2$C1_v_C2_FDR[i] < 0.05 ~ paste0(tb[i,5],"*"),
    .default = paste0(tb[i,5])
  )
}
wb2 <- createWorkbook()
addWorksheet(wb2,sheetName = "Sheet1")
writeDataTable(wb2,sheet = "Sheet1",tb,tableStyle = "none")
addStyle(wb2,"Sheet1",SigCliffStyle,cols = 2,rows = 1+which(tobold[,1]))
addStyle(wb2,"Sheet1",SigCliffStyle,cols = 3,rows = 1+which(tobold[,2]))
addStyle(wb2,"Sheet1",SigCliffStyle,cols = 4,rows = 1+which(tobold[,3]))
addStyle(wb2,"Sheet1",SigCliffStyle,cols = 5,rows = 1+which(tobold[,4]))
saveWorkbook(wb2,"Clean_res.xlsx",overwrite = T)

##- (4) Vizualizations
plotData <- ResultsTable2
plotData <- plotData[!is.na(plotData$Sig_description),]

ggplot(plotData, aes(Sig_description,fill = Subsystem))+
  geom_bar(position = "dodge")+
  coord_flip()+
  theme_classic()
ggsave("SigTypes.svg",width = 24,height = 12)


plotData <- ResultsTable2 |> 
  select(c(Task_ID,
           Task_Name,
           Control_v_DCM_FDR,Control_v_DCM_FC,
           Control_v_C1_FDR,Control_v_C1_FC,
           Control_v_C2_FDR,Control_v_C2_FC,
           C1_v_C2_FDR,
           C1_v_C2_FC)) |> 
  pivot_longer(cols = c(Control_v_DCM_FDR,Control_v_C1_FDR,Control_v_C2_FDR,C1_v_C2_FDR),
               names_to = "Contrast1",values_to = "FDR") |> 
  pivot_longer(cols = c(Control_v_DCM_FC,Control_v_C1_FC,Control_v_C2_FC,C1_v_C2_FC),
               names_to = "Contrast2",values_to = "FC") |> 
  mutate(Contrast1 = str_remove(Contrast1,"_FDR")) |> 
  mutate(Contrast2 = str_remove(Contrast2,"_FC")) |>
  filter(Contrast1 == Contrast2) |> 
  mutate(Significance = FDR < 0.05 & (FC >= 1.1 | FC <= -1.1)) |> 
  filter(Significance) |> 
  mutate(Task_Name = str_trim(Task_Name)) |> 
  mutate(Task_Name = str_replace(Task_Name,"( .*? .*? .*?) ","\\1\n")) |> 
  mutate(Task_Name = str_remove(Task_Name,"\n$")) |> 
  mutate(TOI = Task_ID %in% c("Task_347","Task_272","Task_321","Task_18",
                              "Task_322","Task_271","Task_351","Task_79",
                              "Task_350","Task_273","Task_324","Task_352")) |> 
  mutate(Task_Name = ifelse(TOI,glue("<b>{Task_Name}<b>"),Task_Name))
# 
# library(scales)
# squish_trans <- function(from, to, factor) {
#   
#   trans <- function(x) {
#     
#     if (any(is.na(x))) return(x)
#     
#     # get indices for the relevant regions
#     isq <- x > from & x < to
#     ito <- x >= to
#     
#     # apply transformation
#     x[isq] <- from + (x[isq] - from)/factor
#     x[ito] <- from + (to - from)/factor + (x[ito] - to)
#     
#     return(x)
#   }
#   
#   inv <- function(x) {
#     
#     if (any(is.na(x))) return(x)
#     
#     # get indices for the relevant regions
#     isq <- x > from & x < from + (to - from)/factor
#     ito <- x >= from + (to - from)/factor
#     
#     # apply transformation
#     x[isq] <- from + (x[isq] - from) * factor
#     x[ito] <- to + (x[ito] - (from + (to - from)/factor))
#     
#     return(x)
#   }
#   
#   # return the transformation
#   return(trans_new("squished", trans, inv))
# }

textsize <- 14
g1 <- plotData |> 
  filter(Contrast1 == "Control_v_DCM") |> 
  mutate(AbsFC = abs(FC)) |> 
  arrange(-AbsFC) |> 
  slice_head(n = 10) |> 
  mutate(Ord = dense_rank(FC)) |>
  arrange(Ord) |> 
  mutate(Task_Name = factor(Task_Name, levels = Task_Name[Ord])) |> 
  ggplot(aes(FC,Task_Name)) +
  geom_point() +
  ggbreak::scale_x_break(c(-1.15,1.3),space = 0.2,ticklabels=c(-1.7,-1.5, 1.5)) +
  scale_x_continuous(breaks = c(-1.7,-1.5, 1.5,1.7)) +
  
  theme_pubr(base_size = textsize) +
 # theme(axis.text.y = element_text(face = ifelse(plotData$TOI,"bold","plain")))+
  ylab("Task name") +
  ggtitle("DCM vs Controls")

g2 <- plotData |> 
  filter(Contrast1 == "Control_v_C1") |> 
  mutate(AbsFC = abs(FC)) |> 
  arrange(-AbsFC) |> 
  slice_head(n = 10) |> 
  mutate(Ord = dense_rank(FC)) |>
  arrange(Ord) |> 
  mutate(Task_Name = factor(Task_Name, levels = Task_Name[Ord])) |> 
  ggplot(aes(FC,Task_Name)) +
  geom_point() +
  ggbreak::scale_x_break(c(-1.15,1.15),space = 0.2) +
  scale_x_continuous(breaks = c(-1.4,-1.2, 1.2,1.5)) +
  theme_pubr(base_size = textsize) +
    ggtitle("Metabotype 1 vs Controls") +
  ylab("Task name")

g3 <- plotData |> 
  filter(Contrast1 == "Control_v_C2") |> 
  mutate(AbsFC = abs(FC)) |> 
  arrange(-AbsFC) |> 
  slice_head(n = 10) |> 
  mutate(Ord = dense_rank(FC)) |>
  arrange(Ord) |> 
  mutate(Task_Name = factor(Task_Name, levels = Task_Name[Ord])) |> 
  ggplot(aes(FC,Task_Name)) +
  geom_point() +
  ggbreak::scale_x_break(c(-1.15,1.15),space = 0.2) +
  scale_x_continuous(breaks = c(-2.2,-1.5,-1.2, 1.2,1.5)) +
  theme_pubr(base_size = textsize) + 
  ggtitle("Metabotype 2 vs Controls") +
  ylab("Task name")
  
g4 <- plotData |> 
  filter(Contrast1 == "C1_v_C2") |> 
  mutate(AbsFC = abs(FC)) |> 
  arrange(-AbsFC) |> 
  slice_head(n = 10) |> 
  mutate(Ord = dense_rank(FC)) |>
  arrange(Ord) |> 
  mutate(Task_Name = factor(Task_Name, levels = Task_Name[Ord])) |> 
  ggplot(aes(FC,Task_Name)) +
  geom_point() +
  ggbreak::scale_x_break(c(-1.15,1.15),space = 0.2) +
  scale_x_continuous(breaks = c(-2.2,-1.5,-1.2, 1.2,1.5)) +
  theme_pubr(base_size = textsize) + 
  ggtitle("Metabotype 2 vs Metabotype 1") +
  theme(axis.text.y = ggtext::element_markdown()) +
  ylab("Task name")
g4

g5 <- (g1 | g2)/(g3 | g4)

# TaskName <- plotData |> 
#   filter(Contrast1 == "Control_v_DCM") |> 
#   mutate(AbsFC = abs(FC)) |> 
#   arrange(-AbsFC) |> 
#   slice_head(n = 15)
# plotData |> 
#   mutate(AbsFC = abs(FC)) |> 
#   arrange(-AbsFC) |> 
#   slice_head(n = 15) |> 
# g5 <- g1


ggsave("Plots/topFC_clusts.png",g5,dpi = 300,width = 32,height = 30, units = "cm")
#############"

save.image("../../Data/saved_wkspaces/6_end.RData")
#-- Initial examination revealed that tasks only sig in one cluster have abysmal effect size even for that comparison
