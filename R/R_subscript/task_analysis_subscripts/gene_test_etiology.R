gene_test_etiology <- function(gene,geneName,gxData,phenoData,toSave = T){
  dat <- log(unlist(gxData[gxData$gene == gene,-1]))
  plotData <- data.frame(Gene_expr = dat,
                         Group = phenoData$etiology)
  plotData$Group[plotData$Group == "NF"] <- "Control"
  ggplot(plotData,aes(x=Group,y=Gene_expr))+
    geom_boxplot(aes(fill = Group))
  ttestDat <- t.test(dat[phenoData$etiology == "NF"],dat[phenoData$etiology == "DCM"])
  wilcox.test(dat[phenoData$etiology == "NF"],dat[phenoData$etiology == "DCM"])
  Dat2 <- data.frame(gene = dat,
                     Etiology = phenoData$etiology,
                     Sex = phenoData$gender,
                     Ethnicity = phenoData$race,
                     Age = phenoData$age)
  summary(lm(gene ~ .,Dat2))
  library(ggprism)
  library(glue)
  pvalDat <- data.frame(group1 = "Control",
                        group2 = "DCM",
                        p.adj = ttestDat$p.value,
                        y.position = max(plotData$Gene_expr)*1.05)
  p <- ggplot(plotData,aes(x=Group,y=Gene_expr))+
    geom_boxplot(aes(fill = Group))+
    ylab("log(geTMM) value")+
    ggtitle(paste0(geneName," expression"))+
    theme_cardio()+
    add_pvalue(pvalDat,label = "P-value = {format(p.adj,digits = 3,scientific = -2)}")
  print(p)
  if(toSave){
    ggsave(paste0(geneName,"_plot_etiology.svg"),p,width = 12,height = 12,units = "cm")
  }
  return(p)
}
