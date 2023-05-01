#Check enrichemnt of differantially expressed genes in WGCNA modules
#Scripts were derived from https://github.com/konopkalab/Hippo_Subfields

library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)
enableWGCNAThreads()

#load GeneSet and name of the heatmap
wgcna = list.files(pattern = 'ModuleOutput*')
tab=read.table(wgcna, sep="\t",header=T)
colnames(tab)=c("Gene","DEFINITION","Kwithin")
Genes=as.data.frame(table(tab$DEFINITION))
module_enrich <- "Module_Enrich_Fisher_DEG.RData"
enrich_heatmap <- "Heatmap_GeneSets_Fisher_coolmod_adj_DEG.png"

#load GeneSet and name of the heatmap for Inhibitory
load("DG_GeneSet_InN.RData")
a <- c("InN19_UP", "InN19_DOWN")
GeneSets <- GeneSet
heatmap.label = c("Upregulated", "Downregulated")

##############################################################################################
##################################Enrichment analysis#########################################
# The loop goes along the length of the GeneSets lists
for(i in 1:length(GeneSets)){
  GeneSets[[i]] <- GeneSets[[i]][GeneSets[[i]]$Gene %in% tab$Gene,]
}
ln=length(GeneSets) # number of geneset
cl=length(Genes$Var1) # number of module
md = as.data.frame(as.character(Genes$Var1)) # list of module
md$Freq = as.numeric(rep(0,cl)) # list of module add o as initial start
colnames(md) = c("Module", "Freq")
md <- as.data.frame(md)

TEMP=list()
INT=list()
for (i in 1:ln){
  TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
  INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
  test <- as.data.frame(INT[[i]])
  tmp <- merge(md, test, by.x="Module", by.y="Var1", sort = FALSE, all.x = TRUE)
  INT[[i]]= tmp %>% rowwise() %>% mutate("Freq"= sum(c_across(Freq.x:Freq.y))) %>% mutate_all(~replace(., is.na(.), 0))  %>% dplyr::select(c(Module, Freq)) %>% arrange(Module)
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)
#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT)){
  INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
  INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq
  INT[[i]]$d <- 34160-Genes$Freq-INT[[i]]$Freq
}

# sum(Genes$Freq)
RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value),
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
}
# run
FisherMat=list()
for (i in 1:length(INT)){
  FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
  rownames(FisherMat[[i]]) <- INT[[i]]$Module
  FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)
save(FisherMat,file= module_enrich)

tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT)){
  tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
  FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT)){
  tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,6]))
  FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Pval Adjusted
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=length(a)) 
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames
FisherAdj[FisherAdj>0.05]=1
FisherOR[FisherOR < 1]=0
# annotation for DEG
colnames(FisherAdj)= a
colnames(FisherOR)= a

#get coolmodule list
traitP = list.files(pattern = 'modTraitP_*')
trait=read.table(traitP,sep="\t",header=T)
coolmod <- trait %>% filter(Condition.bi < 0.05) %>% row.names() %>% as.data.frame()
colnames(coolmod) <-"V1"
coolmod <- sub("ME", "", coolmod$V1)

#if all the module
coolmod <- md$Module
FisherPf=FisherAdj[rownames(FisherAdj) %in% coolmod,]
FisherORf=FisherOR[rownames(FisherOR) %in% coolmod,]
df=-log10(FisherPf)
LabelMat = paste(signif(FisherORf, 2), "\n(",signif(FisherPf, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
rownames(df) <- paste("ME", rownames(df), sep="")
png(enrich_heatmap, width = 400, height= 850)
par(mar = c(8, 9, 3, 3))
labeledHeatmap(Matrix = df, xLabels = heatmap.label,  yLabels = rownames(df), colorLabels =TRUE, ySymbols = rownames(df), colors=colorRampPalette(c("white", "#60AF8B"))(50), textMatrix=LabelMat, setStdMargins = FALSE, cex.text = 0.7, xLabelsAngle = 45, main = paste("DG Modules GeneSets Enrichment"), cex.lab = 1)
dev.off()
enrich_heatmap
