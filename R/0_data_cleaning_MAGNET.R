#-------------------------------------------------#
#   0- Cleaning data from the magnet consortium   #
#-------------------------------------------------#

#Libraries
library(edgeR)
library(biomaRt)
library(ggplot2)
#Loading data
rawCounts <- read.csv(file = "../Data/magnet_raw_counts.csv")
colnames(rawCounts)[1] <- "GeneName"
phenoData <- read.csv(file = "../Data/magnet_phenoData.csv")
geneLength <- read.table(file = "../Data/magnet_GeneLengths.txt", sep = "\t",header = TRUE)

#-----------#
# Cleaning  #
#-----------#

#Removing unwanted samples
phenoData <- phenoData[phenoData$etiology == "DCM"|phenoData$etiology == "NF",]
#phenoData <- phenoData[phenoData$gender == "Male",]
#phenoData <- phenoData[phenoData$race == "Caucasian",]
#phenoData <- phenoData[phenoData$Diabetes == "No" | is.na(phenoData$Diabetes),]

colnames(phenoData)[1] <- "ID"
phenoData <- phenoData[-c(2,9,10,11,12,16,18,19,20)]
phenoData$Hypertension[is.na(phenoData$Hypertension)] <- "No"
phenoData$Diabetes[is.na(phenoData$Diabetes)] <- "No"

#Removing patients who fit the legal definition of dwarfism
phenoData <- phenoData[phenoData$height > 142 | is.na(phenoData$height),]
#calculating BMI
phenoData$BMI <- phenoData$weight/(phenoData$height/100)^2


phenoDataDCM <- phenoData[phenoData$etiology == "DCM",-2]
phenoDataControl <- phenoData[phenoData$etiology == "NF",-2]
write.csv(phenoDataDCM,"../Results/0_Clean_Data/magnet_phenoData_DCM.csv",row.names = FALSE)

rawCounts <- rawCounts[,c(TRUE,colnames(rawCounts)[-1] %in% phenoData$ID)]
rownames(rawCounts) <- rawCounts$GeneName
rawCounts <- rawCounts[geneLength$ENSEMBL,]

geneLength$geneLengthK <- geneLength$genelength /1000
#Performing rpk normalization
rpkCounts <- rawCounts

for (i in 2:ncol(rpkCounts)){
  rpkCounts[,i] <- rpkCounts[,i]/geneLength$geneLengthK
}

### -geTMM for all samples- ###

geTMMCounts <- DGEList(rpkCounts[,-1])
geTMMCounts <- edgeR::calcNormFactors(geTMMCounts,method = "TMM")
geTMMCounts <- edgeR::cpm(geTMMCounts)
geTMMCounts <- log(geTMMCounts + 1)
plot(colSums(geTMMCounts))
#boxplot(geTMMCounts,outline = FALSE)

pcaRes <- pcaMethods::pca(t(geTMMCounts))


###May be removed, removing sample with lowest expression before re-normalizing
isMin <- which.min(colSums(geTMMCounts))
isMax <- which.max(colSums(geTMMCounts))
extreme <- rep("No",nrow(pcaRes@scores))
extreme[isMin] <- "Min"
extreme[isMax] <- "Max"
plotData <- data.frame(extreme = extreme,
                       PC1 = pcaRes@scores[,1],
                       PC2 = pcaRes@scores[,2],
                       libSize = colSums(geTMMCounts),
                       condition = phenoData$etiology,
                       batch = phenoData$Library.Pool,
                       batch_cond = paste0(phenoData$Library.Pool,"_",phenoData$etiology))
ggplot(plotData,aes(PC1,PC2, color = condition))+
  geom_point()+
  ggtitle("All")
ggplot(plotData,aes(PC1,PC2, color = extreme))+
  geom_point()+
  ggtitle("All")
ggplot(plotData,aes(PC1,PC2, color = libSize))+
  geom_point()+
  ggtitle("All")
ggplot(plotData,aes(PC1,PC2, color = batch_cond))+
  geom_point()+
  ggtitle("All")


geTMM_BC <- removeBatchEffect(geTMMCounts,phenoData$Library.Pool)


pcaRes_BC <- pcaMethods::pca(t(geTMM_BC))

plotData_BC <- data.frame(extreme = extreme,
                       PC1 = pcaRes_BC@scores[,1],
                       PC2 = pcaRes_BC@scores[,2],
                       libSize = colSums(geTMM_BC),
                       condition = phenoData$etiology,
                       batch = phenoData$Library.Pool,
                       batch_cond = paste0(phenoData$Library.Pool,"_",phenoData$etiology))
ggplot(plotData_BC,aes(PC1,PC2, color = condition))+
  geom_point()+
  ggtitle("All_BC")
ggplot(plotData_BC,aes(PC1,PC2, color = extreme))+
  geom_point()+
  ggtitle("All_BC")
ggplot(plotData_BC,aes(PC1,PC2, color = libSize))+
  geom_point()+
  ggtitle("All_BC")
ggplot(plotData_BC,aes(PC1,PC2, color = batch_cond))+
  geom_point()+
  ggtitle("All_BC")

# After examination, no samples to remove
geTMM_BC_unlog <- exp(geTMM_BC)
sum(geTMM_BC_unlog < 0)
toSave <- data.frame(gene = rownames(geTMM_BC),as.data.frame(geTMM_BC_unlog))
#colnames(toSave)[1] <- "gene"
#Remove min also in phenoData
phenoData <- phenoData[phenoData$ID %in% colnames(toSave),]

libPool <- phenoData$Library.Pool
phenoData[,11] <- NULL
phenoData$LibPool <- libPool

write.csv(toSave,"../Results/0_Clean_Data/geTMM_MAGNET_BC.csv",row.names = FALSE)
write.csv(phenoData,"../Results/0_Clean_Data/magnet_phenoData_complete.csv",row.names = FALSE)

### -Normalizing for DCM samples- ###
rpkCounts <- rawCounts
rpkCounts  <- rpkCounts[,c(T,phenoData$etiology == "DCM")]
for (i in 2:ncol(rpkCounts)){
  rpkCounts[,i] <- rpkCounts[,i]/geneLength$geneLengthK
}
rpkCounts <- rpkCounts
geTMMCounts <- DGEList(rpkCounts[,-1])
geTMMCounts <- edgeR::calcNormFactors(geTMMCounts,method = "TMM")
geTMMCounts <- edgeR::cpm(geTMMCounts)
# 
# geTMMCounts <- DGEList(rpkCounts[,-1])
# geTMMCounts <- edgeR::calcNormFactors(geTMMCounts,method = "TMM")
# geTMMCounts <- edgeR::cpm(geTMMCounts)


geTMMCounts <- log(geTMMCounts+1)
plot(colSums(geTMMCounts+1))
#boxplot(geTMMCounts,outline = FALSE)

pcaRes <- pcaMethods::pca(t(geTMMCounts))
isMin <- which.min(colSums(geTMMCounts))
isMax <- which.max(colSums(geTMMCounts))
extreme <- rep("No",nrow(pcaRes@scores))
extreme[isMin] <- "Min"
extreme[isMax] <- "Max"
plotData <- data.frame(extreme = extreme,
                       PC1 = pcaRes@scores[,1],
                       PC2 = pcaRes@scores[,2],
                       libSize = colSums(geTMMCounts),
                       batch = phenoDataDCM$Library.Pool)

ggplot(plotData,aes(PC1,PC2, color = extreme))+
  geom_point()+
  ggtitle("DCM")
ggplot(plotData,aes(PC1,PC2, color = libSize))+
  geom_point()+
  ggtitle("DCM")


#Batch effect removal
geTMM_DCM_BC <- removeBatchEffect(geTMMCounts,phenoDataDCM$Library.Pool)

pcaRes <- pcaMethods::pca(t(geTMM_DCM_BC))
isMin <- which.min(colSums(geTMM_DCM_BC))
isMax <- which.max(colSums(geTMM_DCM_BC))
extreme <- rep("No",nrow(pcaRes@scores))
extreme[isMin] <- "Min"
extreme[isMax] <- "Max"
plotData <- data.frame(extreme = extreme,
                       PC1 = pcaRes@scores[,1],
                       PC2 = pcaRes@scores[,2],
                       libSize = colSums(geTMM_DCM_BC),
                       batch = phenoDataDCM$Library.Pool)

ggplot(plotData,aes(PC1,PC2, color = extreme))+
  geom_point()+
  ggtitle("DCM_BC")
ggplot(plotData,aes(PC1,PC2, color = libSize))+
  geom_point()+
  ggtitle("DCM_BC")
ggplot(plotData,aes(PC1,PC2, color = batch))+
  geom_point()+
  ggtitle("DCM_BC")


toSave <- data.frame(rownames(geTMM_DCM_BC),as.data.frame(exp(geTMM_DCM_BC)))
colnames(toSave)[1] <- "gene"
#boxplot(unlist(geTMMCounts))
write.csv(toSave,"../Results/0_Clean_Data/geTMM_MAGNET_BC_DCM.csv",row.names = FALSE)
libPool <- phenoDataDCM$Library.Pool
phenoDataDCM[,10] <- NULL
phenoDataDCM$LibPool <- libPool
write.csv(phenoDataDCM,"../Results/0_Clean_Data/magnet_phenoData_DCM.csv",row.names = FALSE)


### -Normalizing for Control samples- ###
rpkCounts <- rawCounts
rpkCounts  <- rpkCounts[,c(T,phenoData$etiology == "NF")]
for (i in 2:ncol(rpkCounts)){
  rpkCounts[,i] <- rpkCounts[,i]/geneLength$geneLengthK
}
rpkCounts <- rpkCounts
geTMMCounts <- DGEList(rpkCounts[,-1])
geTMMCounts <- edgeR::calcNormFactors(geTMMCounts,method = "TMM")
geTMMCounts <- edgeR::cpm(geTMMCounts)
# 
# geTMMCounts <- DGEList(rpkCounts[,-1])
# geTMMCounts <- edgeR::calcNormFactors(geTMMCounts,method = "TMM")
# geTMMCounts <- edgeR::cpm(geTMMCounts)


geTMMCounts <- log(geTMMCounts+1)
plot(colSums(geTMMCounts+1))
#boxplot(geTMMCounts,outline = FALSE)

pcaRes <- pcaMethods::pca(t(geTMMCounts))
isMin <- which.min(colSums(geTMMCounts))
isMax <- which.max(colSums(geTMMCounts))
extreme <- rep("No",nrow(pcaRes@scores))
extreme[isMin] <- "Min"
extreme[isMax] <- "Max"
plotData <- data.frame(extreme = extreme,
                       PC1 = pcaRes@scores[,1],
                       PC2 = pcaRes@scores[,2],
                       libSize = colSums(geTMMCounts),
                       batch = phenoDataControl$Library.Pool,
                       sex = phenoDataControl$gender)

ggplot(plotData,aes(PC1,PC2, color = extreme))+
  geom_point()+
  ggtitle("Control")
ggplot(plotData,aes(PC1,PC2, color = libSize))+
  geom_point()+
  ggtitle("Control")
ggplot(plotData,aes(PC1,PC2, color = batch))+
  geom_point()+
  ggtitle("Control")
ggplot(plotData,aes(PC1,PC2, color = sex))+
  geom_point()+
  ggtitle("Control")


#Batch effect removal
geTMM_NF_BC <- removeBatchEffect(geTMMCounts,phenoDataControl$Library.Pool)

pcaRes <- pcaMethods::pca(t(geTMM_NF_BC))
isMin <- which.min(colSums(geTMM_NF_BC))
isMax <- which.max(colSums(geTMM_NF_BC))
extreme <- rep("No",nrow(pcaRes@scores))
extreme[isMin] <- "Min"
extreme[isMax] <- "Max"
plotData <- data.frame(extreme = extreme,
                       PC1 = pcaRes@scores[,1],
                       PC2 = pcaRes@scores[,2],
                       libSize = colSums(geTMM_NF_BC),
                       batch = phenoDataControl$Library.Pool,
                       sex = phenoDataControl$gender,
                       age = phenoDataControl$age,
                       diabetes = phenoDataControl$Diabetes,
                       BMI = phenoDataControl$weight / (phenoDataControl$height/100)^2,
                       hypertension = phenoDataControl$Hypertension,
                       LVEF = phenoDataControl$LVEF,
                       ethnicity = phenoDataControl$race)

ggplot(plotData,aes(PC1,PC2, color = extreme))+
  geom_point()+
  ggtitle("Control_BC")
ggplot(plotData,aes(PC1,PC2, color = libSize))+
  geom_point()+
  ggtitle("Control_BC")
ggplot(plotData,aes(PC1,PC2, color = batch))+
  geom_point()+
  ggtitle("Control_BC")
ggplot(plotData,aes(PC1,PC2, color = ethnicity))+
  geom_point()+
  ggtitle("Control_BC")
ggplot(plotData,aes(PC1,PC2, color = sex))+
  geom_point()+
  ggtitle("Control_BC")


toSave <- data.frame(rownames(geTMM_NF_BC),as.data.frame(exp(geTMM_NF_BC)))
colnames(toSave)[1] <- "gene"
#boxplot(unlist(geTMMCounts))
write.csv(toSave,"../Results/0_Clean_Data/geTMM_MAGNET_BC_NF.csv",row.names = FALSE)
write.csv(phenoDataControl[,-10],"../Results/0_Clean_Data/magnet_phenoData_NF.csv",row.names = FALSE)


