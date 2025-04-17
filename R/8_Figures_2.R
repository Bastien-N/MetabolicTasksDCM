#....................................#
#             Libraries              #
#....................................#

library(tidyverse)
library(ggpubr)
# Supplementary materials, histogram of shapiro-test p-value ----
#....................................#
##           Loading data             #----
#....................................#
score_all <- read.table("../Results/2_MAGNET_NFvDCM/Cellfie_scores/cellfie_score.txt",
                        sep = ",")


sample_table <- read.table("../Results/0_Clean_Data/magnet_phenoData_complete.csv",sep = ",",header = T)
task_table <- read.table("../Results/0_Clean_Data/task_report.txt",sep = ",",header = F)
isMinusOne <- apply(score_all,1,\(x){any(x == -1)})
#....................................#
##           Reformatting             #----
#....................................#
task_table <- task_table[!task_table$V1 %in% c("277","258"),]
task_table <- task_table[-c(351,352),]
task_table <- task_table[!isMinusOne,]
score_all <- score_all[!isMinusOne,]





#....................................#
##             Analysis               #----
#....................................#
score_all |> 
  t() |> 
  as.data.frame() |> 
  filter(sample_table$etiology != "NF") |> 
  map(\(x) {
    Shapiro_p = shapiro.test(x)$p.value
  }) |> 
  unlist() |> 
  data.frame(p_value = _) |> 
  mutate(log_p_value = log(p_value)) |> 
  ggplot(aes(-log_p_value)) +
  geom_histogram() +
  geom_vline(xintercept = -log(0.05))
ggsave("../Results/2_MAGNET_NFvDCM/Plots/Shapiro_pval.png",dpi = "retina")