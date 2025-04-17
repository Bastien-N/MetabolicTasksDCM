# Libraries
library(tidyverse)
library(ggpubr)
library(ggprism)
library(pheatmap)
library(stringr)
# Data loading
subjTable <- read.table("../Results/0_Clean_Data/magnet_phenoData_complete.csv", sep = ",", header = T)
clustTable <- read.table("../Results/3_MAGNET_DCM_clustering/Tables/DCM_clusts_from_DCM.txt",header = T)

subjTable$Cluster <- clustTable$Cluster[match(subjTable$ID,clustTable$Sample)]
subjTable$Cluster[is.na(subjTable$Cluster)] <- "Control"
subjTable$Cluster <- str_replace(subjTable$Cluster,"_"," ")
subjTable$Cluster <- factor(subjTable$Cluster, levels = c("Control","Cluster 1","Cluster 2"))

mwRes_all <- read.table("../Results/4_MAGNET_clusts_v_NF/Tables/mann_whitney_tasks_clusters_v_control.txt", header = T)
mwRes_clust <- read.table("../Results/3_MAGNET_DCM_clustering/Tables/mann_whitney_cluster_tasks.txt", header = T)
scores_all <- read.table("../Results/2_MAGNET_NFvDCM/Cellfie_scores/cellfie_score.txt",sep = ",") |> 
  filter(V1 != -1) |> 
  magrittr::set_colnames(subjTable$ID) |> 
  add_column(Task = unique(mwRes_all$Task_ID),.before = 1)
taskTable <- read.table("../Results/0_Clean_Data/task_report.txt",sep = ",") |> 
  mutate(Task_ID = paste0("Task_",V1))
taskTable <- taskTable[!taskTable$Task_ID %in% c("Task_277","Task_258"),]
taskTable <- taskTable[taskTable$Task_ID %in% mwRes_all$Task_ID,]




# Making the TOI plot ----
isTOI <- taskTable$Task_ID %in% c("Task_347","Task_272","Task_321","Task_18",
                                  "Task_322","Task_271","Task_351","Task_79",
                                  "Task_350","Task_273","Task_324","Task_352")
g1 <- scores_all[isTOI,] |> 
  as.data.frame(row.names = 1:sum(isTOI)) |> 
  column_to_rownames("Task") |> 
  t() |> 
  as.data.frame()  |> 
  mutate(Group = subjTable$Cluster) |> 
  pivot_longer(cols = 1:12,names_to = "Task",values_to = "Score")  |> 
  mutate(TaskName = taskTable$V4[match(Task,taskTable$Task_ID)]) |> 
  mutate(TaskName = str_replace_all(TaskName,"-phosphate","-p")) |> 
  mutate(TaskName = str_remove(TaskName,"\\(.*\\)")) |> 
  mutate(TaskName = ifelse(str_count(TaskName," ") > 2,
                           paste0(str_extract(TaskName,"[^ ]* [^ ]* [^ ]* "),"\n",
                                  str_remove(TaskName,"[^ ]* [^ ]* [^ ]* ")),
                           TaskName)) |> 
  ggplot(aes(Group,Score))+
  #geom_violin() +
   # stat_compare_means(comparisons = comps,method = "wilcox.test",
   #                    label = "p.signif",hide.ns = T,vjust = 0.5) +
  facet_wrap(vars(TaskName),scales = "free_y") +
  
  geom_point(size = 0.5,position = position_dodge2(width = 0.5),
             shape = 1,fill = NA,alpha = 0.8) +
  geom_boxplot(fill = NA,outlier.shape = NA) +
  theme_pubr() +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 35,vjust = 0.1)) +
  ylab("Task activity score")
g1
ggsave("../Results/4_MAGNET_clusts_v_NF/Plots/TOI_figure.png",g1,
       height = 21,width = 21,units = "cm",
       dpi = "retina")