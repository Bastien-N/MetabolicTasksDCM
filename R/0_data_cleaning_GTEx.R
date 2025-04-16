#--------------------------------#
#   0- Cleaning data from GTEx   #
#--------------------------------#


library(CePa)
library(data.table)
library(ineq)
library(stringr)
#library(refGenome)
library(edgeR)
library(pcaMethods)
library(ggplot2)
library(glmnet)

gtexData <- fread("../Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct",skip = 2)
genome <- read.table("../data/gencode.v26.GRCh38.genes.gtf",skip = 6,header = FALSE,sep = "\t")
sampleData <- fread("../Data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
subjectData <- fread("../Data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
#examining samples
table(sampleData$SMTSD)

#excluding samples without overlap in gtexData and sampleData4
samples1 <- colnames(gtexData)
samples1 <- samples1[-c(1,2)]

samples2 <- sampleData$SAMPID

samples2 <- samples2[samples2 %in% samples1]
sampleData <- sampleData[sampleData$SAMPID %in% samples2,]
gtexData <- gtexData[,-2]
gtexData <- as.data.frame(gtexData)
#Reducing to samples of interest
groups <- as.factor(sampleData$SMTSD)

organsOfInterest <- c("Adipose - Subcutaneous","Colon - Transverse","Heart - Left Ventricle","Liver","Muscle - Skeletal")
samplesOI <- sampleData$SMTSD %in% organsOfInterest
sampleData <- sampleData[samplesOI,]

gtexData <- gtexData[,c(TRUE,samplesOI)]



#Processing genome file
genome <- genome[genome$V3 == "gene",]

colnames(genome) <- c("Chr","Source","Feature","start","end","score","strand","frame","group")
genome_gene_id <- str_extract(genome$group,"ENSG[:digit:]*\\.[:digit:]*[_PARY]*")
genome$transcriptID <- genome_gene_id
genome$geneID <- str_remove(genome$transcriptID,"\\..*")

dataTranscripts <- gtexData$Name
toKeep <- !str_detect(dataTranscripts,"PAR_")
dataTranscripts <- dataTranscripts[toKeep]

dataGenes <- str_remove(dataTranscripts,"\\..*")
gtexData <- gtexData[toKeep,-1]

#Getting genes that are in the model
modelGenes <- unlist(read.table("../Data/model_genes.txt",header = F))
isInModel <- dataGenes %in%modelGenes



rpkCounts <- gtexData
genome <- genome[genome$transcriptID %in% dataTranscripts,]
rownames(genome) <- genome$transcriptID
genome <- genome[dataTranscripts,]

geneLengths <- (genome$end - genome$start) /1000

rpkCounts <- sweep(rpkCounts,1,geneLengths,`/`)

#Examining possible filtering

sum(rowSums(rpkCounts ==0) == ncol(rpkCounts))
sum(rowSums(rpkCounts ==0 ) == ncol(rpkCounts) & !isInModel)


#Not filtering to avoid issues

geTMMCounts <- DGEList(rpkCounts)
geTMMCounts <- edgeR::calcNormFactors(geTMMCounts,method = "TMM")
geTMMCounts <- edgeR::cpm(geTMMCounts)


geTMMCounts <- log(geTMMCounts+1)
hist(geTMMCounts)
pcaRes <- pcaMethods::pca(t(geTMMCounts))

plotData <- data.frame(PC1 = pcaRes@scores[,1],
                       PC2 = pcaRes@scores[,2],
                       Tissue = sampleData$SMTSD,
                       Batch = sampleData$SMGEBTCH,libSize = colSums(geTMMCounts),
                       Center = sampleData$SMCENTER)


ggplot(plotData,aes(PC1,PC2,color = Tissue))+
  geom_point()
ggsave("gtex_before_correction.png")
ggplot(plotData,aes(PC1,PC2,color = libSize))+
  geom_point()
ggsave("gtex_before_correction_libSize.png")
ggplot(plotData,aes(PC1,PC2,color = Center))+
  geom_point()
ggsave("gtex_before_correction_center.png")
ggplot(plotData,aes(PC1,PC2,color = Batch))+
  geom_point()+
  theme(legend.position = "none")
ggsave("gtex_before_correction_batches.png")
toSave <- data.frame(gene = dataGenes,
                     exp(geTMMCounts)-1)


write.csv(toSave,"../Results/0_Clean_Data/geTMM_GTEx_noBC.csv",row.names = FALSE)

geTMMCounts_BC <- removeBatchEffect(geTMMCounts,plotData$Batch) 

pcaRes <- pcaMethods::pca(t(geTMMCounts_BC))

plotData <- data.frame(PC1 = pcaRes@scores[,1],
                       PC2 = pcaRes@scores[,2],
                       Tissue = sampleData$SMTSD,
                       Batch = sampleData$SMGEBTCH,libSize = colSums(geTMMCounts_BC),
                       Center = sampleData$SMCENTER)


ggplot(plotData,aes(PC1,PC2,color = Tissue))+
  geom_point()
ggsave("gtex_after_correction.png")
ggplot(plotData,aes(PC1,PC2,color = libSize))+
  geom_point()
ggsave("gtex_after_correction_libSize.png")
ggplot(plotData,aes(PC1,PC2,color = Center))+
  geom_point()
ggsave("gtex_after_correction_center.png")
ggplot(plotData,aes(PC1,PC2,color = Batch))+
  geom_point()+
  theme(legend.position = "none")
ggsave("gtex_after_correction_batches.png")

toSave <- data.frame(gene = dataGenes,
                     exp(geTMMCounts_BC))
write.csv(toSave,"../Results/0_Clean_Data/geTMM_GTEx_BC.csv",row.names = FALSE)


#examining batches distribution
chitest <- chisq.test(table(sampleData$SMGEBTCH,sampleData$SMTSD))
chitest$p.value #batches are not homogeneous among tissue types

#linear regression
makeDummy <- function(x){
  groups <- unique(x)
  res <- matrix(nrow = length(x),ncol = length(groups))
  colnames(res) <- groups
  for(i in 1:length(groups)){
    for (ii in 1:length(x)){
      if(x[ii] == groups[i]){
        res[ii,i] <- 1
      }else{
        res[ii,i] <- 0
      }
    }
  }
  return(res)
}
accuracyByClass <- function(pred,real){
  classes <- unique(real)
  res <- data.frame(class = classes, percent_correct = rep(0,length(classes)))
  for(i in 1:length(classes)){
    res[i,2] <- sum(pred == real & real == classes[i])  / sum(real == classes[i])
  }
  return(res)
}
batchDummy <- makeDummy(sampleData$SMGEBTCH)
tissueDummy <- makeDummy(sampleData$SMTSD)
linearDat <- data.frame(tissue = tissueDummy[,1],
                        batchDummy)
logModel <- glmnet::glmnet(batchDummy,sampleData$SMTSD,
                           family = "multinomial",
                           alpha = 1)

predTest <- predict(logModel,batchDummy,type = "class")

accuracy <- apply(predTest,2,function(p){
  sum(sampleData$SMTSD == p)
})
table(sampleData$SMTSD,sampleData$SMTSD == predTest[,40])
accuracyByClass(predTest[,40],sampleData$SMTSD)

############################
## Same but using only cardiac data
gxCardiac <- rpkCounts[,sampleData$SMTSD == "Heart - Left Ventricle"]
sampleCardiac <- sampleData[sampleData$SMTSD == "Heart - Left Ventricle",]

geTMMCardiac <- DGEList(gxCardiac)
geTMMCardiac <- edgeR::calcNormFactors(geTMMCardiac,method = "RLE")
geTMMCardiac <- edgeR::cpm(geTMMCardiac)


geTMMCardiac <- log(as.matrix(geTMMCardiac)+1)
hist(geTMMCardiac)
pcaRes <- pcaMethods::pca(t(geTMMCardiac))


plotData <- data.frame(PC1 = pcaRes@scores[,1],
                       PC2 = pcaRes@scores[,2],
                       Tissue = sampleCardiac$SMTSD,
                       Batch = sampleCardiac$SMGEBTCH,libSize = colSums(geTMMCardiac),
                       Center = sampleCardiac$SMCENTER,
                       method = sampleCardiac$SMNABTCHT)


# ggplot(plotData,aes(PC1,PC2,color = Tissue))+
#   geom_point()
#ggsave("gtex_before_correction.png")
ggplot(plotData,aes(PC1,PC2,color = libSize))+
  geom_point()
ggplot(plotData,aes(PC1,PC2,color = Center))+
  geom_point()
ggplot(plotData,aes(PC1,PC2,color = method))+
  geom_point()

#ggsave("gtex_before_correction_libSize.png")
ggplot(plotData,aes(PC1,PC2,color = Batch))+
  geom_point()+
  theme(legend.position = "none")
#ggsave("gtex_before_correction_batches.png")

geTMMCardiac_BC <- removeBatchEffect(geTMMCardiac,plotData$Center)

pcaRes <- pcaMethods::pca(t(geTMMCardiac_BC))

plotData <- data.frame(PC1 = pcaRes@scores[,1],
                       PC2 = pcaRes@scores[,2],
                       Tissue = sampleCardiac$SMTSD,
                       Batch = sampleCardiac$SMGEBTCH,libSize = colSums(geTMMCardiac_BC),
                       Center = sampleCardiac$SMCENTER)


ggplot(plotData,aes(PC1,PC2,color = Tissue))+
  geom_point()
#ggsave("gtex_before_correction.png")
ggplot(plotData,aes(PC1,PC2,color = libSize))+
  geom_point()
ggplot(plotData,aes(PC1,PC2,color = Center))+
  geom_point()

#ggsave("gtex_before_correction_libSize.png")
ggplot(plotData,aes(PC1,PC2,color = Batch))+
  geom_point()+
  theme(legend.position = "none")

###- phenoData construction -###

subject <- sampleData$SAMPID
subject <- str_remove(subject,"-.{4}-SM-.*$")
subjectIndex <- rep(0,nrow(sampleData))
for(i in 1:length(subjectIndex)){
  subjectIndex[i] <- which(str_detect(subjectData$SUBJID,subject[i]))
}
cod <- as.character(subjectData$DTHHRDY[subjectIndex])
cod[is.na(cod)] <- "5"
cod <- sapply(cod,function(co){
  switch(co,
         "0" = "Ventilator",
         "1" = "Violent fast",
         "2" = "Natural fast",
         "3" = "Intermediate",
         "4" = "Slow",
         "5" = "Unknown")})

sex <- sapply(as.character(subjectData$SEX[subjectIndex]),function(s){
  switch(s,
         "1" = "M",
         "2" = "F")
})
phenoData <- data.frame(sampleID = sampleData$SAMPID,
                        subjectID = subject,
                        tissue = sampleData$SMTS,
                        sex = sex,
                        age = as.numeric(str_remove(subjectData$AGE[subjectIndex],"-.{2}$")),
                        cause_death = cod)
write.csv(phenoData,"../Results/0_Clean_Data/GTEx_phenoData.csv",
          row.names = FALSE)
