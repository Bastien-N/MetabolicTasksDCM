gene_test_clusters <- function(gene,geneName,gxdata,group,toSave = T){
  dat3 <- log(unlist(gxData[gxData$gene == gene,-1]))
  
  utest1 <- wilcox.test(dat3[group == "Control"],dat3[group == "Cluster_1"])
  utest2 <- wilcox.test(dat3[group == "Control"],dat3[group == "Cluster_2"])
  utest3 <- wilcox.test(dat3[group == "Cluster_1"],dat3[group == "Cluster_2"])
  pvals <- c(utest1$p.value,utest2$p.value,utest3$p.value)
  qvals <- p.adjust(pvals,"BH")
  
  plotData <- data.frame(Gene_expr = dat3,
                         Group = group)
  pvalDat <- data.frame(group1 = c("Control","Control","Cluster_1"),
                        group2 = c("Cluster_1","Cluster_2","Cluster_2"),
                        p.adj = qvals,
                        label = paste0("FDR = ",format(qvals,digits = 3,scientific = -2)),
                        y.position = max(plotData$Gene_expr)*c(1.0,1.05,1.1))
  isSig <- pvalDat$p.adj < 0.05
  p <- ggplot(plotData,aes(x=Group,y=Gene_expr))+
    geom_boxplot(aes(fill = Group))+
    geom_violin(fill = NA,linewidth = 0.1)+
    ylab("log(geTMM) value")+
    ggtitle(paste0(geneName," expression"))+
    theme_cardio()
  
    if(sum(!isSig) == length(isSig)){
      p <- p + add_pvalue(pvalDat[!isSig,],tip.length = 0.001,linetype = "dashed" )
    }else if(sum(!isSig) == 0){
      p <- p + add_pvalue(pvalDat[isSig,],tip.length = 0.001,fontface = "bold")
    }else{
      p <- p + add_pvalue(pvalDat[!isSig,],tip.length = 0.001,linetype = "dashed" )+
      add_pvalue(pvalDat[isSig,],tip.length = 0.001,fontface = "bold")
    }
    
  print(p)
  if(toSave){
    ggsave(paste0(geneName,"_plot_clusters.svg"),p,width = 12,height = 12,units = "cm")
  }
  return(p)
}
