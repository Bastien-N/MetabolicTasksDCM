#------------------#
#   TITLE          #
#   BSC Nihant     #
#------------------#

#....................................#
#             Libraries              #----
#....................................#

library(tidyverse)
library(CIBERSORT)
library(biomaRt)
library(rstatix)
library(ggpubr)
library(patchwork)
#....................................#
#           Loading data             #----
#....................................#
rna_data <- read.table("../Data/magnet_raw_counts.csv",sep = ",",header = T)
geneLength <- read.table("../Data/magnet_GeneLengths.txt",header = T)
reference <- CIBERSORT::LM22

groups <- read.table("../Results/4_MAGNET_DCM_clustering/Tables/DCM_clusts_from_DCM.txt",
                     header  =T)
sampleData <- read.table("../Results/0_Clean_Data/magnet_phenoData_complete.csv",
                         header = T,
                         sep = ",")
#....................................#
#           Reformatting             #----
#....................................#
sampleData$Group <- sampleData$etiology
for (i in 1:nrow(sampleData)){
  if (sampleData$etiology[i] == "DCM") {
    sampleData$Group[i] <- groups$Cluster[groups$Sample == sampleData$ID[i]]
  }
}
geneLength <- geneLength[rowSums(rna_data[,-1]) > 1,]
rna_data <- rna_data[rowSums(rna_data[,-1]) > 1 ,c("X",sampleData$ID)]
mart <- useMart("ensembl")
mart <- useMart("ensembl","hsapiens_gene_ensembl")
# eMart <- useEnsembl("ensembl","hsapiens_gene_ensembl",
#                     mirror = "useast")

biomart_res <- getBM(c("ensembl_gene_id","external_gene_name"),
                     filters = "external_gene_name",
                     values = rownames(reference),
                     mart = mart)

biomart_res2 <- getBM(c("ensembl_gene_id","external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = rna_data$X,
                     mart = mart)

fpkm_data <- rna_data
scaling <- rowSums(fpkm_data[,2:ncol(fpkm_data)]) * geneLength[,2]
dat <- fpkm_data[,2:ncol(fpkm_data)]
fpkm_data[,2:ncol(fpkm_data)] <- 1e9 * dat / scaling
rownames(fpkm_data) <- fpkm_data[,1]
fpkm_data <- fpkm_data[,-1]
fpkm_data <- fpkm_data[!is.na(fpkm_data[,2]),]

fpkm_symbol <- fpkm_data |> 
  mutate(Symbol = biomart_res2$external_gene_name[match(row.names(fpkm_data),biomart_res2$ensembl_gene_id)]) |> 
  filter(!is.na(Symbol) & Symbol != "") |> 
  group_by(Symbol) |> 
  summarise_all(sum)
fpkm_symbol <- as.data.frame(fpkm_symbol)
row.names(fpkm_symbol) <- fpkm_symbol$Symbol
fpkm_symbol <- as.matrix(fpkm_symbol[,-1])


ref_ensembl <- reference[match(biomart_res$external_gene_name,rownames(reference)),]
rownames(ref_ensembl) <- biomart_res$ensembl_gene_id 

fpkm_symbol |> 
  as.data.frame() |> 
  rownames_to_column("Gene") |> 
  relocate(Gene,.before = 1) |> 
  filter(Gene %in% row.names(reference)) |> 
write.table(x = _, "../Results/0_Clean_Data/magnet_mixture.txt",
            row.names = T,
            sep = "\t")

reference |> 
  as.data.frame() |> 
  rownames_to_column("GeneSymbol") |> 
  relocate(GeneSymbol,.before = 1) |> 
write.table(x = _, "../Results/0_Clean_Data/LM22.csv",
            row.names = T,
            sep = "\t")

#....................................#
#             Analysis               #----
#....................................#
sampleData$Group[sampleData$Group == "NF"] <- "Control"
sampleData$Group <- factor(sampleData$Group,
                           levels = c("Control","Cluster_1","Cluster_2"))




# Ensembl

test_ciber <- cibersort(sig_matrix = ref_ensembl,
                        mixture_file = as.matrix(fpkm_data),
                        perm = 10000,
                        QN = F)

# test_ciber3 <- cibersort(sig_matrix = ref_ensembl,
#                         mixture_file = as.matrix(fpkm_data),
#                         perm = 0)
plotData <- as.data.frame(test_ciber)
plotData$Groups <- factor(sampleData$Group,
                          levels = c("Control","Cluster_1","Cluster_2"))

comps <- list(c("Control","Cluster_1"),
              c("Control","Cluster_2"),
              c("Cluster_1","Cluster_2"))
for (i in 1:22){
  g1 <- ggplot(plotData,aes(Groups,!!sym(colnames(plotData)[i])))+
    geom_violin() +
    stat_compare_means(comparisons = comps,method = "wilcox.test",
                       label = "p.format") +
    geom_boxplot(fill = NA) +
    theme_pubr()
  
  chiTest <- chisq_test(table(plotData$Groups, 0 != plotData[,i]))
  
  g2 <- ggplot(plotData,aes(Groups, fill = 0 != !!sym(colnames(plotData)[i])))+
    geom_bar(position = position_dodge()) +
    labs(fill = paste0(colnames(plotData)[i]," present")) +
    annotate("text",label = paste0("p.val = ",format(chiTest$p,scientific = -1)),
             x = 1.2,y = 75) +
    theme_pubr()
    
  ggsave(paste0("../Results/11_deconvolution/Plots/Attemp_1/ensembl_",colnames(plotData)[i],".png"),
         g1 + g2,width = 21,height = 12,units = "cm")
}

# Symbol
test_ciber2 <- cibersort(sig_matrix = reference,
                         mixture_file = fpkm_symbol,
                         perm = 10000,
                         QN = F)
plotData <- as.data.frame(test_ciber2)
plotData$Groups <- factor(sampleData$Group,
                          levels = c("Control","Cluster_1","Cluster_2"))

comps <- list(c("Control","Cluster_1"),
              c("Control","Cluster_2"),
              c("Cluster_1","Cluster_2"))
for (i in 1:22){
  g1 <- ggplot(plotData,aes(Groups,!!sym(colnames(plotData)[i])))+
    geom_violin() +
    stat_compare_means(comparisons = comps,method = "wilcox.test",
                       label = "p.format") +
    geom_boxplot(fill = NA) +
    theme_pubr()
  
  chiTest <- chisq_test(table(plotData$Groups, 0 != plotData[,i]))
  
  g2 <- ggplot(plotData,aes(Groups, fill = 0 != !!sym(colnames(plotData)[i])))+
    geom_bar(position = position_dodge()) +
    labs(fill = paste0(colnames(plotData)[i]," present")) +
    annotate("text",label = paste0("p.val = ",format(chiTest$p,scientific = -1)),
             x = 1.2,y = 75) +
    theme_pubr()
  
  ggsave(paste0("../Results/11_deconvolution/Plots/Attemp_2/symbol_",colnames(plotData)[i],".png"),
         g1 + g2,width = 21,height = 12,units = "cm")
}

save.image("../Data/saved_wkspaces/11_Deconvolution.RData")


# Loading CibersortX results
test_ciberx <- read.table("../Results/11_deconvolution/CIBERSORTx_Job1_Results.txt",
                          header = T,
                          sep = "\t")

plotData <- test_ciberx
plotData$Groups <- factor(str_replace(sampleData$Group,"_"," "),
                          levels = c("Control","Cluster 1","Cluster 2"))

comps <- list(c("Control","Cluster 1"),
              c("Control","Cluster 2"),
              c("Cluster 1","Cluster 2"))
for (i in 2:23){
  g1 <- ggplot(plotData,aes(Groups,!!sym(colnames(plotData)[i])))+
    geom_violin() +
    stat_compare_means(comparisons = comps,method = "wilcox.test",
                       label = "p.format") +
    geom_boxplot(fill = NA) +
    theme_pubr()
  
  chiTest <- chisq_test(table(plotData$Groups, 0 != plotData[,i]))
  
  g2 <- ggplot(plotData,aes(Groups, fill = 0 != !!sym(colnames(plotData)[i])))+
    geom_bar(position = position_dodge()) +
    labs(fill = paste0(colnames(plotData)[i]," present")) +
    annotate("text",label = paste0("p.val = ",format(chiTest$p,scientific = -1)),
             x = 1.2,y = 75) +
    theme_pubr()
  
  ggsave(paste0("../Results/11_deconvolution/Plots/Attemp_x/symbol_",colnames(plotData)[i],".png"),
         g1 + g2,width = 21,height = 12,units = "cm")
}

g1 <- plotData |> 
  pivot_longer(cols = c(B.cells.naive:Neutrophils),
               names_to = "Celltype",
               values_to = "EstimatedProp") |> 
  mutate(Celltype = str_replace_all(Celltype,"\\."," ")) |> 
  mutate(Celltype = str_replace(Celltype,"T cell","T-cell")) |> 
  mutate(Celltype = str_replace(Celltype,"B cell","B-cell")) |>
  mutate(Celltype = str_remove(Celltype,"  Tregs")) |> 
  mutate(Celltype = ifelse(str_count(Celltype," ") > 1,
         paste0(str_extract(Celltype,".* .* "),"\n",
                str_remove(Celltype,".* .* ")),
         Celltype)) |> 
  ggplot(aes(Groups,EstimatedProp))+
  #geom_violin() +
  stat_compare_means(comparisons = comps,method = "wilcox.test",
                     label = "p.signif",hide.ns = T,vjust = 0.5) +
  facet_wrap(vars(Celltype),scales = "free_y") +
  
  geom_point(size = 0.5,position = position_dodge2(width = 0.5),
             shape = 1,fill = NA,alpha = 0.8) +
  geom_boxplot(fill = NA,outlier.shape = NA) +
  theme_pubr() +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 35,vjust = 0.1)) +
  ylab("Estimated cell-type proportion")
g1
ggsave("../Results/11_deconvolution/Plots/Attemp_x/all_celltypes.png",g1,
       height = 21,width = 19,units = "cm",
       dpi = "retina")
