library(Seurat)

harmonized_object <- readRDS(file = 
	"/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

Idents(harmonized_object) <- harmonized_object$RNA_snn_res.0.7

#Add BRETIGEA module scores

print("Adding BRETIGEA scores")
library(BRETIGEA)

markers_for_annot <- split(markers_df_human_brain$markers, markers_df_human_brain$cell)


#Assess numbers of overlapping genes 

lapply(markers_for_annot, function(x){
	c(length(x), length(intersect(x, rownames(harmonized_object))))
})

print("Adding BRETIGEA module scores")
harmonized_object <- AddModuleScore(harmonized_object, markers_for_annot, nbin = 10, ctrl = 500, search = TRUE)
for(i in 1:length(markers_for_annot)){
	item_name <- paste("Cluster", i, sep ="")
	harmonized_object@meta.data[names(markers_for_annot)[i]] <- harmonized_object@meta.data[item_name]
	harmonized_object@meta.data[item_name] <- NULL
}

#Do same thing as above with UCell to AddModuleScores

library(UCell)

harmonized_object <- AddModuleScore_UCell(harmonized_object, features = markers_for_annot, chunk.size = 500)

#plot BRETIGEA

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/2_BRETIGEA_scores_Seurat_UCell.pdf", height = 16, width =12)
VlnPlot(harmonized_object, features = names(markers_for_annot), pt.size = -1, ncol=2)
VlnPlot(harmonized_object, features = paste0(names(markers_for_annot), "_UCell"), pt.size = -1, ncol = 2)
dev.off()

#Add BrainInABlender module scores

library(BrainInABlender)

print("Adding BrainInABlender scores")

#This includes mouse genes 
BIAB_markers <- split(CellTypeSpecificGenes_Master3$GeneSymbol_Human, CellTypeSpecificGenes_Master3$CellType_Primary)
BIAB_markers <- lapply(BIAB_markers, na.omit)
BIAB_markers <- lapply(BIAB_markers, as.character)


lapply(BIAB_markers, function(x){
        c(length(x),length(intersect(x, rownames(harmonized_object))))
})

print("Adding BrainInABlender module scores")
harmonized_object <- AddModuleScore(harmonized_object, BIAB_markers, nbin = 10)
for(i in 1:length(BIAB_markers)){
	item_name <- paste("Cluster", i, sep = "")
	harmonized_object@meta.data[names(BIAB_markers)[i]] <- harmonized_object@meta.data[item_name]
	harmonized_object@meta.data[item_name] <- NULL
}

#Do same thing as above with UCell

harmonized_object <- AddModuleScore_UCell(harmonized_object, features = BIAB_markers, chunk.size = 500)

#plot BrainInABlender

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/2_BrainInABlender_scores_Seurat_UCell.pdf", height = 16, width =14)
VlnPlot(harmonized_object, features = names(BIAB_markers), pt.size = -1, ncol=2)
VlnPlot(harmonized_object, features = paste0(names(BIAB_markers), "_UCell"), pt.size = -1, ncol=2)
dev.off()

#write out the module scores
library(readr)

write_csv(harmonized_object@meta.data[,c(names(markers_for_annot), names(BIAB_markers), paste0(names(BIAB_markers), "_UCell"), paste0(names(markers_for_annot), "_UCell"))], 
	file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/2_BRETIGEA_BrainInABlender_scores.csv")

#find markers Seurat

#apparently this step below will take 12 hours, could try with future if required 

#all_markers <- FindAllMarkers(harmonized_object, test.use = "poisson", latent.vars = c("Batch", "Sample", "Chemistry", "percent.mt", "nCount_RNA"))
#write_csv(all_markers, file = "/home/malosree/def-gturecki/projects/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/all_markers_Seurat_poisson.csv")

#find markers Presto

print("Finding presto markers")

library(tidyverse)

all_markers_presto <- presto::wilcoxauc(harmonized_object, "RNA_snn_res.0.7")
all_markers_presto <- filter(as_tibble(all_markers_presto), padj < 0.05 & logFC > log(1.5) & pct_in-pct_out > 10)

write_csv(all_markers_presto, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/2_all_markers_Presto_WilcoxAUC.csv")

#Run GOs

print("Running GOs")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

all_markers_presto_split <- split(all_markers_presto, all_markers_presto$group)
egos <- list()
onts <- c("CC", "BP", "MF")
for(ont in onts) {
	egos[[ont]] <- lapply(all_markers_presto_split, function(x){
		print(length(x$feature))
		enrichGO(gene = x$feature, universe= rownames(harmonized_object), OrgDb = org.Hs.eg.db, ont = ont, pAdjustMethod="BH", 
			pvalueCutoff=0.01, qvalueCutoff=0.05, keyType = 'SYMBOL')
	})
}

for(ont in onts) {
	print("creating plots")
	plots_list <- lapply(names(egos[[ont]]), function(x){
		this_clust <- egos[[ont]][[x]]
		if(length(this_clust$Description) > 0) barplot(this_clust)+ggtitle(paste("Cluster", x))
		else return("No over-representation")
	})
	pdf(file = paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/2_GO_plot_Presto_markers_", 
	ont, ".pdf"),
		


onefile = TRUE, height = 7, width = 10)
	print("saving plots")
	for(res in plots_list) {
		if(typeof(res) == "character") message(res)
		else print(res)
	}	
	dev.off()
}

#make plots

library(ggplot2)

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/2_Unlabelled_DotPlot_main.pdf",
	height = 12, width = 10)
DotPlot(harmonized_object, features = c("SNAP25", "RBFOX3", "SLC17A7", "SATB2", "GAD1", "GAD2", 
		"ALDH1L1", "GFAP", "GJA1", "MBP", "PLP1", "PDGFRA", "VIM", "CX3CR1"))+
		theme(axis.text.x = element_text(angle = 90))
dev.off()

#Layer specific marker gene expression

print("Layer specificity assessment")

#Load layer specificity data from Maynard et al.
Maynard_layer_spec <- read_csv(file = "/home/malosree/projects/def-gturecki/malosree/41593_2020_787_MOESM3_ESM_Table_S5.csv")

print("Filter by FDR")
Maynard_layer_spec %>% filter(fdr < 0.01) -> Maynard_layer_spec

#select the clusters with overlapping markers
print("Filter markers data")
all_markers_presto %>% filter(feature %in% Maynard_layer_spec$gene_name) %>% group_by(feature) -> overlap_markers

#select the layer information for the markers detected in our data
print("combine layer info")
Layer_info <- Maynard_layer_spec %>% dplyr::select(gene_name, Layer1, Layer2, Layer3, Layer4, Layer5, Layer6, WM) %>% 
		filter(gene_name %in% overlap_markers$feature) %>% 
		group_by(gene_name) %>%
		summarise_all(sum) -> Layer_info

print("Combine with marker information")
#combine the information
full_join(overlap_markers, Layer_info, by = c("feature" = "gene_name")) %>% 
	ungroup() %>% group_by(group) %>% 
	write_csv(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/2_cluster_Maynard_layer_info.csv")

print(summary(warnings()))
sessionInfo()
