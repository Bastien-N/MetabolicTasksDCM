#---------------------------------------------#
#   6- Examining Actual gene data from tasks  #
#---------------------------------------------#

##- (1) Setup -#
localDir <- "../Results/6_DE_clusters/"
if(!dir.exists(localDir)){ dir.create(localDir) }
setwd(localDir)
if(!dir.exists("Tables")){ dir.create("Tables") }
if(!dir.exists("Plots")){ dir.create("Plots")}
#-- Libraries
library(biomaRt)
library(edgeR)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggpubr)
library(topGO)
library(stringr)
library(dplyr)
library(tidyr)
library(magrittr)
library(rrvgo)
source("../../R/R_subscript/task_analysis_load.R")

#-- Data
gxData0 <- read.table(file = "../../Data/magnet_raw_counts.csv",sep = ",",header = T)
gxData_DCM_TMM <- read.table(file = "../0_Clean_Data/geTMM_MAGNET_BC_DCM.csv",sep = ",",header = T)
clusterTable <- read.table("../3_MAGNET_DCM_clustering/Tables/DCM_clusts_from_DCM.txt",header = T)
phenoData <- read.csv("../0_Clean_Data/magnet_phenoData_complete.csv")


AMPK_genes <- unlist(read.table("../../Data/AMPK_genes.txt"))
vATPase_genes <- unlist(read.table("../../Data/v_ATPase_genes.txt"))


modelGenes <- read.table("../../Data/model_genes.txt",header = F) %>% unlist()
mart <- useEnsembl("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl")

trad <- getBM(attributes = c("ensembl_gene_id",
                             "entrezgene_id",
                             "gene_biotype"),
              filters = "ensembl_gene_id",
              values = gxData0$X,
              mart = mart)
org <- org.Hs.eg.db
keytypes(org)

#-- Functions
prepare_go <- function(DEtest) {
  inGo <- DEtest$table$logFC
  
  names(inGo) <- DEtest$genes$ID
  inGo <- sort(inGo, decreasing = T)
  return(inGo)
}
prepare_testGO <- function(DEres){
  require(topGO)
  require(org.Hs.eg.db)
  inGo <- DEres$table$logFC
  
  names(inGo) <- DEres$genes$ID
  inGo <- sort(inGo, decreasing = T)
  minsig <- max(DEres$table$logFC[DEres$table$logFC < 0 & DEres$table$PValue < 0.05])
  maxsig <- min(DEres$table$logFC[DEres$table$logFC > 0 & DEres$table$PValue < 0.05])
  selection <- function(x) {
    x > maxsig | x < minsig
  }
  allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="ensembl")
  GOdata <- new("topGOdata", ontology="BP", allGenes=inGo,
                annot=annFUN.GO2genes, GO2genes=allGO2genes,
                geneSel=selection, nodeSize=5,
                description = "GO data")
  return(GOdata)
}
gen_GO_results <- function(GOdata,ksRes) {
  goEnrichment <- GenTable(GOdata, KS=ksRes, orderBy="KS", topNodes=100)
  goEnrichment$KS <- as.numeric(goEnrichment$KS)
  goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
  goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
  goEnrichment
}

better_go_plot <- function(goTest,title) {
  g <- dotplot(goTest,
               x = "GeneRatio",
               title = title)
  lev <- levels(g$data$Description)
  wds <- strsplit(lev, "\\W+")
  for(i in 1:length(wds)) {
    if(length(wds[[i]]) >= 6 ){
      lev[i] <- paste0(str_remove(lev[i],paste0(wds[[i]][6],".*$")),"...")
    }
  }
  levels(g$data$Description) <- lev
  ggplot(g$data,aes(GeneRatio,Description, size = Count)) + 
    geom_point() + 
    theme_pubr() +
    theme(panel.grid.major.y = element_line(color = "grey"),
          plot.title = element_text(hjust = 0.5)) +
    ylab("Gene ontology") +
    ggtitle(title)
}

prepare_kegg <- function(DEtest,BM) {
  trad2 <- getBM(attributes = c("ensembl_gene_id",
                                "entrezgene_id",
                                "gene_biotype"),
                 filters = "ensembl_gene_id",
                 values = DEtest$genes$ID,
                 mart = BM)
  
  inKegg <- DEtest$table$logFC
  names(inKegg) <- paste0("TMP",names(inKegg))
  toKeep <- rep(F,length(inKegg))
  eIDs <- rep("",length(inKegg))
  for (i in 1:length(inKegg)) {
    g <- DEtest$genes$ID[i]
    eID <- trad2$entrezgene_id[trad2$ensembl_gene_id == g]
    if( sum(trad2$ensembl_gene_id == g) == 0){next()}
    if(length(eID) > 1){eID <- eID[1]}
    if(!any(is.na(eID)) & length(eID) > 0) {
      toKeep[i] <- T
      vals <- inKegg[!is.na(trad2$entrezgene_id) & trad2$entrezgene_id == eID]
      if (length(vals) > 1) {
        toKeep[!is.na(trad2$entrezgene_id) & trad2$entrezgene_id == eID] <- abs(vals) == max(abs(vals))
      }
    }
    eIDs[i] <- as.character(eID)
  }
  inKegg <- inKegg[toKeep]
  eIDs <- eIDs[toKeep]
  names(inKegg) <- eIDs
  
  inKegg <- sort(inKegg,decreasing = TRUE)
  return(inKegg)
}
##- (2) Building model
colnames(gxData0)[1] <- "gene"
gxData0 <- gxData0[,colnames(gxData0) %in% c("gene",phenoData$ID)]

gxData0 <- gxData0[,c(1,1 + order(colnames(gxData0)[-1]))]


phenoData <- phenoData[order(phenoData$ID),]
group <- sapply(phenoData$ID,\(id){
  if(id %in% clusterTable$Sample) {
    return(as.character(clusterTable$Cluster[clusterTable$Sample == id]))
  } else {
    return("Control")
  }
})
group <- factor(group, levels = c("Control", "Cluster_1", "Cluster_2"))
phenoData$etiology <- group2 <- factor(phenoData$etiology,levels = c("NF","DCM"))

gxData <- edgeR::DGEList(gxData0[,-1],
                         samples = phenoData,
                         group = group,
                         genes = data.frame(ID = gxData0$gene))
gxData <- gxData[filterByExpr(gxData),,keep.lib.size = F]

gxData <- calcNormFactors(gxData, method = "TMM")
#-- For DCM v NF
design <- model.matrix( ~ LibPool + group, data = gxData$samples)
design
gxData <- estimateDisp(gxData,design = design)
fit <- glmQLFit(gxData,design = design)

save.image("../../Data/saved_wkspaces/6_DE_clusters_modelBuilt.RData")
##- (3) DE analysis C1 v C2
test_clust <- glmTreat(fit,contrast = makeContrasts(C1 = "groupCluster_2 - groupCluster_1",levels = design),lfc = log2(1.5))
topTags(test_clust)

inGo_clust <- prepare_go(test_clust)
goTest_clust <- gseGO(inGo_clust,
                      ont = "BP",
                      eps = 0,
                      org,
                      keyType = "ENSEMBL",
                      pvalueCutoff = 0.01)
better_go_plot(goTest_clust, title = "Top 10 GO terms between clusters")

ggsave("Plots/gsea_go_C1_C2.png")
write.table(goTest_clust@result[,1:9],file = "Tables/gsea_go_C1_C2.txt",sep = "\t",row.names = F)



inKegg_clust <- prepare_kegg(test_clust,mart)
keggTest_clust <- gseKEGG(inKegg_clust,eps = 1e-15,pAdjustMethod = "BH",pvalueCutoff = 0.2)
better_go_plot(keggTest_clust, title = "Top 10 Kegg terms between clusters")
ggsave("Plots/gsea_kegg_C1_C2.png")
write.table(keggTest_clust@result[,1:9],file = "Tables/gsea_kegg_C1_C2.txt",sep = "\t",row.names = F)
##- (4) DE analysis Control v C1
test_c1 <- glmTreat(fit,coef = 13,lfc = log2(1.5))
topTags(test_c1)

inGo_c1 <- prepare_go(test_c1)
goTest_c1 <- gseGO(inGo_c1,
                      ont = "BP",
                      eps = 0,
                      org,
                      keyType = "ENSEMBL",
                      pvalueCutoff = 0.01
                   )
better_go_plot(goTest_c1, title = "Top 10 GO terms for cluster 1")
ggsave("Plots/gsea_go_Control_C1.png")
write.table(goTest_c1@result[,1:9],file = "Tables/gsea_go_Control_C1.txt",sep = "\t",row.names = F)


inKegg_c1 <- prepare_kegg(test_c1,mart)
keggTest_c1 <- gseKEGG(inKegg_c1,eps = 1e-15,pAdjustMethod = "BH",pvalueCutoff = 0.2)
better_go_plot(keggTest_c1, title = "Top 10 Kegg terms for cluster 1")
ggsave("Plots/gsea_kegg_Control_C1.png")
write.table(keggTest_c1@result[,1:9],file = "Tables/gsea_kegg_Control_C1.txt",sep = "\t",row.names = F)
##- (3) DE analysis Control v C2
test_c2 <- glmTreat(fit,coef = 14,lfc = log2(1.5))
topTags(test_c2)

inGo_c2 <- prepare_go(test_c2)
goTest_c2 <- gseGO(inGo_c2,
                   ont = "BP",
                   eps = 0,
                   org,
                   keyType = "ENSEMBL")
better_go_plot(goTest_c2, title = "Top 10 GO terms for cluster 2")
ggsave("Plots/gsea_go_Control_C2.png")
write.table(goTest_c2@result[,1:9],file = "Tables/gsea_go_Control_C2.txt",sep = "\t",row.names = F)


inKegg_c2 <- prepare_kegg(test_c2,mart)
keggTest_c2 <- gseKEGG(inKegg_c2,eps = 1e-15,pAdjustMethod = "BH",pvalueCutoff = 0.2)
better_go_plot(keggTest_c2, title = "Top 10 Kegg terms for cluster 2")
ggsave("Plots/gsea_kegg_Control_C2.png")
write.table(keggTest_c2@result[,1:9],file = "Tables/gsea_kegg_Control_C2.txt",sep = "\t",row.names = F)


save.image("../../Data/saved_wkspaces/6_DE_clusters_testRun.RData")


##- (4) specific protein complexes

AMPK_Res <- test_clust$table[test_clust$genes$ID %in% AMPK_genes,] |> 
  rbind(... = _,test_c1$table[test_c1$genes$ID %in% AMPK_genes,]) |> 
  rbind(... = _,test_c2$table[test_c2$genes$ID %in% AMPK_genes,]) |> 
  mutate(Gene = rep(test_clust$genes$ID[test_clust$genes$ID %in% AMPK_genes],3)) |> 
  mutate(Cmp = c(rep("C1_v_C2",length(Gene)/3),
                 rep("Cont_v_C1",length(Gene)/3),
                 rep("Cont_v_C2",length(Gene)/3))) |> 
  extract(,c(5,6,1,2,3,4))
write.table(AMPK_Res,"Tables/DE_AMPK.txt",sep = "\t",row.names = F)

vATPase_Res <- test_clust$table[test_clust$genes$ID %in% vATPase_genes,] |> 
  rbind(... = _,test_c1$table[test_c1$genes$ID %in% vATPase_genes,]) |> 
  rbind(... = _,test_c2$table[test_c2$genes$ID %in% vATPase_genes,]) |> 
  mutate(Gene = rep(test_clust$genes$ID[test_clust$genes$ID %in% vATPase_genes],3)) |> 
  mutate(Cmp = c(rep("C1_v_C2",length(Gene)/3),
                 rep("Cont_v_C1",length(Gene)/3),
                 rep("Cont_v_C2",length(Gene)/3))) |> 
  extract(,c(5,6,1,2,3,4))
write.table(AMPK_Res,"Tables/DE_vATPase.txt",sep = "\t",row.names = F)

##- (5) Specific pathway


##- (6) rrvgo

simMat_clust <- rrvgo::calculateSimMatrix(goTest_clust@result$ID[goTest_clust@result$qvalue <= 0.01],
                                          orgdb = "org.Hs.eg.db",
                                          ont = "BP")
vgo_score_clust <- -log(goTest_clust@result$pvalue[goTest_clust@result$qvalue <= 0.01])
names(vgo_score_clust) <- goTest_clust@result$ID[goTest_clust@result$qvalue <= 0.01]

redSim_clust <- rrvgo::reduceSimMatrix(simMat_clust,
                                       vgo_score_clust,
                                       orgdb = "org.Hs.eg.db")
png("Plots/treemap_C1_v_C2.png",width = 40,height = 21,"cm",res = 300,bg = "white")

rrvgo::treemapPlot(redSim_clust,
                   fontsize.labels = c(40,0),
                   title = "Metabotype 2 vs Metabotype 1",
                   inflate.labels = T,
                   lowerbound.cex.labels = 1)
graphics.off()


##
simMat_c1 <- rrvgo::calculateSimMatrix(goTest_c1@result$ID[goTest_c1@result$qvalue <= 0.05],
                                          orgdb = "org.Hs.eg.db",
                                          ont = "BP")
vgo_score_c1 <- -log(goTest_c1@result$pvalue[goTest_c1@result$qvalue <= 0.05])
names(vgo_score_c1) <- goTest_c1@result$ID[goTest_c1@result$qvalue <= 0.05]

redSim_c1 <- rrvgo::reduceSimMatrix(simMat_c1,
                                    scores = vgo_score_c1,
                                    orgdb = "org.Hs.eg.db")
png("Plots/treemap_Control_v_C1.png",width = 30,height = 15,"cm",res = 300,bg = "white")
rrvgo::treemapPlot(redSim_c1,
                   title = "Metabotype 1 vs Controls",
                   fontsize.labels = c(10,0),
                   inflate.labels = T)

graphics.off()
##
simMat_c2 <- rrvgo::calculateSimMatrix(goTest_c2@result$ID[goTest_c2@result$qvalue <= 0.01],
                                       orgdb = "org.Hs.eg.db",
                                       ont = "BP")
vgo_score_c2 <- -log(goTest_c2@result$pvalue[goTest_c2@result$qvalue <= 0.01])
names(vgo_score_c2) <- goTest_c2@result$ID[goTest_c2@result$qvalue <= 0.01]

redSim_c2 <- rrvgo::reduceSimMatrix(simMat_c2,
                                    scores = vgo_score_c2,
                                    orgdb = "org.Hs.eg.db")

png("Plots/treemap_Control_v_C2.png",width = 30,height = 15,"cm",res = 300,bg = "white")
rrvgo::treemapPlot(redSim_c2,
                   title = "Metabotype 2 vs Controls",
                   fontsize.labels = c(10,0),
                   inflate.labels = T)
graphics.off()
png("Plots/treemap_Control_v_C2_biosb.png",width = 15,height = 7,"cm",res = 300,bg = "white")
rrvgo::treemapPlot(redSim_c2,
                   title = "Metabotype 2 vs Controls",
                   fontsize.labels = c(10,0),
                   inflate.labels = T)
graphics.off()
# Volcano plot
axis_title_size <- 20
text_size <- 18
test_clust$table |> 
  mutate(p.adjusted = p.adjust(PValue,"BH")) |> 
  mutate(Significance = case_when(
    p.adjusted <= 0.05 & logFC > 0 ~ "Up-regulated",
    p.adjusted <= 0.05 & logFC < 0 ~ "Down-regulated",
    p.adjusted > 0.05 ~ "Non-significant"
  )) |> 
ggplot(aes(logFC,1-log2(p.adjusted), color = Significance)) +
  geom_point() +
  scale_colour_manual(values = c("Up-regulated" = "red",
                    "Down-regulated" = "blue",
                    "Non-significant" = "grey")) +
  geom_vline(xintercept = log2(1.5),linetype = 2,color = "grey") +
  geom_vline(xintercept = -log2(1.5),linetype = 2,color = "grey") +
  theme_pubr() +
  theme(
    # plot.background = element_rect(fill = "#ECFCFF",colour = "#ECFCFF"),
    #     panel.background = element_rect(fill = "#ECFCFF",colour = "#ECFCFF"),
    #     legend.background = element_rect(fill = "#ECFCFF",colour = "#ECFCFF"),
        axis.title.y = element_text(size = axis_title_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.text.x = element_text(size = text_size),
        axis.text.y = element_text(size = text_size),
        legend.text = element_text(size = text_size),
        legend.title = element_text(size = text_size))
ggsave("Plots/volcano_clusts.png",width = 25,height = 12,units = "cm",dpi = 300)
test_c1$table |> 
  mutate(p.adjusted = p.adjust(PValue,"BH")) |> 
  mutate(Significance = case_when(
    p.adjusted <= 0.05 & logFC > 0 ~ "Up-regulated",
    p.adjusted <= 0.05 & logFC < 0 ~ "Down-regulated",
    p.adjusted > 0.05 ~ "Non-significant"
  )) |> 
  ggplot(aes(logFC,1-log2(p.adjusted), color = Significance)) +
  geom_point() +
  scale_colour_manual(values = c("Up-regulated" = "red",
                                 "Down-regulated" = "blue",
                                 "Non-significant" = "grey")) +
  geom_vline(xintercept = log2(1.5),linetype = 2,color = "grey") +
  geom_vline(xintercept = -log2(1.5),linetype = 2,color = "grey") +
  theme_pubr() +
  theme(
    # plot.background = element_rect(fill = "#ECFCFF",colour = "#ECFCFF"),
    #     panel.background = element_rect(fill = "#ECFCFF",colour = "#ECFCFF"),
    #     legend.background = element_rect(fill = "#ECFCFF",colour = "#ECFCFF"),
    axis.title.y = element_text(size = axis_title_size),
    axis.title.x = element_text(size = axis_title_size),
    axis.text.x = element_text(size = text_size),
    axis.text.y = element_text(size = text_size),
    legend.text = element_text(size = text_size),
    legend.title = element_text(size = text_size))
ggsave("Plots/volcano_c1.png",width = 25,height = 12,units = "cm",dpi = 300)

test_c2$table |> 
  mutate(p.adjusted = p.adjust(PValue,"BH")) |> 
  mutate(Significance = case_when(
    p.adjusted <= 0.05 & logFC > 0 ~ "Up-regulated",
    p.adjusted <= 0.05 & logFC < 0 ~ "Down-regulated",
    p.adjusted > 0.05 ~ "Non-significant"
  )) |> 
  ggplot(aes(logFC,1-log2(p.adjusted), color = Significance)) +
  geom_point() +
  scale_colour_manual(values = c("Up-regulated" = "red",
                                 "Down-regulated" = "blue",
                                 "Non-significant" = "grey")) +
  geom_vline(xintercept = log2(1.5),linetype = 2,color = "grey") +
  geom_vline(xintercept = -log2(1.5),linetype = 2,color = "grey") +
  theme_pubr() +
  theme(
    # plot.background = element_rect(fill = "#ECFCFF",colour = "#ECFCFF"),
    #     panel.background = element_rect(fill = "#ECFCFF",colour = "#ECFCFF"),
    #     legend.background = element_rect(fill = "#ECFCFF",colour = "#ECFCFF"),
    axis.title.y = element_text(size = axis_title_size),
    axis.title.x = element_text(size = axis_title_size),
    axis.text.x = element_text(size = text_size),
    axis.text.y = element_text(size = text_size),
    legend.text = element_text(size = text_size),
    legend.title = element_text(size = text_size))
ggsave("Plots/volcano_c2.png",width = 25,height = 12,units = "cm",dpi = 300)

save.image("../../Data/saved_wkspaces/6_end.RData")
#load("../../Data/saved_wkspaces/6_end.RData")



