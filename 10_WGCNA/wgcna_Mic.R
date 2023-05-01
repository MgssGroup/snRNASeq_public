#Conduct WGCNA(https://cran.r-project.org/web/packages/WGCNA/index.html) in Female microglia or Parvalbumin internueon 
#https://github.com/konopkalab/Hippo_Subfields were used as reference

############################
library(muscat)
library(edgeR)
library(limma)
library(tidyverse)
library(WGCNA)
library(tidyr)
library(DESeq2)
library(cluster);
library(reshape2)
library(RColorBrewer)
library(flashClust)
library(SingleCellExperiment)
options(stringsAsFactors = FALSE);

#set files
rds.file <- "Micro_counts.Rds"
meta.file <- "Male_female_metadata_combined.csv"

#get cpm counts and remove subject with low reads count 
res <- readRDS(rds.file)
counts <- as.data.frame(res)
counts<- subset(counts, select=-c(F18)) 
nrow(counts)

#filter for low counts 
filter <- rowSums(counts) > 5
counts_select <- counts[filter,]
nrow(counts_select)

#set covariates in meta file
meta <- read.csv(meta.file, header =T, sep=",")
meta$Condition <- as.factor(meta$Condition)
meta$Batch <- as.character(meta$Batch)
meta$Age <- as.numeric(meta$Age)
meta$pH <- as.numeric(meta$pH)
meta$PMI <- as.numeric(meta$PMI)
datTraits <- as.data.frame(meta)
rownames(datTraits) <- datTraits[,1]

#check format for meta and counts
datTraits <- datTraits %>% filter(Sex=="Female")
datTraits <- datTraits[order(rownames(datTraits)),]  
datTraits <- datTraits[colnames(counts_select),]  
counts_select <- counts_select[,order(colnames(counts_select))]  
all(rownames(datTraits) == colnames(counts_select))

#normalize and transform DESeqDataSetFromMatrix without defining a specific model / variance stabilizing transformation
counts_select <- round(counts_select) %>% as.data.frame()
dds <- DESeqDataSetFromMatrix(countData=counts_select, colData=datTraits, design = ~ 1)
rld <- rlog(dds, blind = TRUE)

#remove batch effect using limma::removeBatcheffect
datTraits$pH <- datTraits$pH %>% replace_na(mean(datTraits$pH, na.rm = T)) #place mean(ph) for the subject with no ph information
cov.matrix <- cbind(datTraits$Age, datTraits$pH, datTraits$PMI)
design0 <- model.matrix(~ Condition, data = datTraits)
normalized <- removeBatchEffect(assay(rld), batch=datTraits$Batch, covariates=cov.matrix , design = design0)
datExpr <- t(normalized)

######################################################################################################################
############################################ Power analysis check ####################################################
powers = c(seq(2,30,2))
sft=pickSoftThreshold(datExpr, powerVector=powers, verbose=5, blockSize= 40000, networkType = "signed")
par(mfrow = c(1,2), mar=c(5.1,5.1,4.1,2.1));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
abline(h=0.5,col="red"); abline(h=0.8,col="blue");abline(h=0.9,col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

######################################################################################################################
############################################ MODULE construction #####################################################
PWR=12
minModule=30
cor <- WGCNA::cor

#check scale free 
adjacency = adjacency(datExpr, power = PWR, type="signed")
TOM = TOMsimilarity(adjacency,TOMType = "signed")
dissTOM = 1-TOM
k=as.vector(apply(adjacency, 2, sum, na.rm=T))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

#build modules
net = blockwiseModules(datExpr, corType="pearson", maxBlockSize = 40000, networkType="signed",power=PWR, minModuleSize= minModule, nThreads=30, TOMType = "signed", TOMDenom = "min", deepSplit=2, verbose=5, mergeCutHeight=0.15, reassignThreshold = 1e-6, detectCutHeight = 0.995, numericLabels=TRUE, saveTOMs=FALSE, pamStage=TRUE, pamRespectsDendro=FALSE)
moduleLabelsAutomatic=net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
MEsAutomatic=net$MEs
unique(moduleColorsAutomatic)
table(moduleColorsAutomatic)
write.table(moduleColorsAutomatic, "DG_colors.txt",sep="\t",quote=F)
save.image(file = "DG_NetData.RData")

#plot dendogram
svg("Dendrogram.Module.svg")
plotDendroAndColors(net$dendrograms[[1]], moduleColorsAutomatic, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

#list hub gene for each module
chooseTopHubInEachModule(datExpr, moduleColorsAutomatic, omitColors = "grey", power = PWR, type = "signed")

#KMEs
KMEs<-signedKME(datExpr, net$MEs,corFnc = "cor", corOptions = "use = 'p'")
kme=data.frame(rownames(counts_select), moduleColorsAutomatic, KMEs)
colnames(kme)[1]="Symbol"
rownames(kme)=NULL
write.table(kme,"KME_DG.txt",sep="\t",quote=F)

#correlate module and phenotype
moduleColorsIEGG=moduleColorsAutomatic
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,moduleColorsIEGG)$eigengenes
MEsIEGG = MEs0
MEsIEGG$MEgrey=NULL

#set phenotype of interest
datTraits$Condition.bi <- gsub("Control",0, datTraits$Condition)
datTraits$Condition.bi <- gsub("Case",1, datTraits$Condition.bi)
datCovariate <- datTraits %>% dplyr::select(Condition.bi, Age, pH, PMI)

modTraitCor= cor(MEsIEGG,datCovariate,method="pearson")
write.table(modTraitCor,"modTraitCor_DG.txt",sep="\t",quote=F)
modTraitP = corPvalueStudent(modTraitCor, nSamples)
write.table(modTraitP,"modTraitP_DG.txt",sep="\t",quote=F)
textMatrix = paste(signif(modTraitCor, 2), "\n(",signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
svg("Heatmap_DatTraits.svg", height=1000, width=700)
par(mar = c(4, 9, 3, 3))
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datCovariate), yLabels = names(MEsIEGG), ySymbols = names(MEsIEGG), colorLabels =FALSE,colors=blueWhiteRed(50),textMatrix=textMatrix, setStdMargins = FALSE, cex.text = 0.6, zlim = c(-1,1), main = paste("Module Association"))
dev.off()

#eigenGeneNet
MEList=moduleEigengenes(datExpr,colors=moduleColorsAutomatic,softPower = PWR,impute = TRUE)
MEs = MEList$eigengenes
MET=orderMEs(MEs)

#GGplotInput
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,moduleColorsAutomatic,softPower = PWR,impute = TRUE)$eigengenes
MEs0$Rows=colnames(counts_select)
MEs0$Class=paste("Class",datTraits$Diagnosis,sep="")
write.table(MEs0, "DG_Matrix_module_correlation.txt",sep="\t",quote=F)

#adjacency matrix
Adj = adjacency(datExpr, power = PWR,type="signed",corFnc = "cor", corOptions = "use = 'p'")
moduleOutput <- data.frame(rownames(counts_select))
moduleOutput[,2]<- moduleColorsAutomatic
intraCon <- intramodularConnectivity(Adj, moduleColorsAutomatic)
moduleOutput[,3]<-intraCon$kWithin
colnames(moduleOutput) <- c("Gene", "ModuleColor", "kWithin")
write.table(moduleOutput, "ModuleOutput_DG.txt", sep="\t", quote=F)
