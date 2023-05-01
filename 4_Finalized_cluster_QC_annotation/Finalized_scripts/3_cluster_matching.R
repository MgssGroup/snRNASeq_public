#label transfer from harmonized male and female data to - male data,  Azimuth, STAB
#Rerun on Dec 7, 2021
#Rerun Sep 15, 2022
#Rerun Apr 17, 2023

#MetaNeighbor 

library(MetaNeighbor)
library(Seurat)
library(SingleCellExperiment) 

#functions for MN comparisons

do_trained_MNUS <- function(dataset_name, study_id, clusters) {
	data <- get(dataset_name)
	study_id <- study_id
	clusters <- clusters
	
	table_file <- paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/MetaNeighbor/3_MNUS_harmonized_to_",
		dataset_name,
		"_",
                clusters)
	
	plot_file <- paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/MetaNeighbor/3_MNUS_harmonized_to_",
		dataset_name,
		"_",
		clusters,
		"_heatmaps")
			
	harmonized_object_clusters <- read.table(paste0(
	"/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/MetaNeighbor/",
	"3_pretrained_harmonized_object_RNA_snn_res.0.7.txt"),
	check.names = FALSE)

	#all hits
	aurocs = MetaNeighborUS(
	trained_model = harmonized_object_clusters, dat = data,
	study_id = data@colData[study_id][[1]], 
	cell_type = data@colData[clusters][[1]],
	fast_version = TRUE
	)

	#best hits
	best_hits = MetaNeighborUS(trained_model = harmonized_object_clusters,
	dat = data,
	study_id = data@colData[study_id][[1]],
	cell_type = data@colData[clusters][[1]],
	one_vs_best = TRUE,
	fast_version = TRUE)

	write.csv(best_hits, file = paste0(plot_file, "_best_hits.csv"))

	#plots
	pdf(file = paste0(plot_file, ".pdf"))
	plotHeatmapPretrained(aurocs, cex = 0.7)
	plotHeatmapPretrained(best_hits, cex = 0.7)
	dev.off()

	svg(file = paste0(plot_file, ".svg"))
	plotHeatmapPretrained(best_hits, cex = 0.7)
	dev.off()
	
	#find and save top hits
	top_hits = topHitsByStudy(aurocs, 
		threshold = 0)

	write.table(top_hits, paste0(table_file, "top_hits_study.csv"))

	#find and save top hits transposed
        top_hits = topHitsByStudy(t(aurocs),
                threshold = 0)

        write.table(top_hits, paste0(table_file, "top_hits_study_transpose.csv"))

} 

#load prep data
harmonized_object <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")
 
#few other misc processing and assessments
 
table(harmonized_object$RNA_snn_res.0.7)
 
#cluster_tree
harmonized_object <- BuildClusterTree(harmonized_object, reduction = "Harmony_Seeded_Batch_Sample_Chemistry", reorder = TRUE)
 
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/3_cluster_tree_ordered.pdf")
PlotClusterTree(harmonized_object)
dev.off()
 
#save the new order of factors, overwriting file from step 1 in this series of analyses
harmonized_object <- StashIdent(harmonized_object, save.name = "RNA_snn_res.0.7")
#commented out on Sep 15, 2022
#saveRDS(harmonized_object, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")
 
#Make a violin plot for Inhib markers and sex-specific genes
library(ggplot2)
 
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/3_inhib_Vlns.pdf",
         height = 12, width = 24)
VlnPlot(harmonized_object, features = c("VIP", "SST", "PVALB", "NPY", "CALB1", "CALB2", "GAD1", "GAD2", "LHX6", "ADARB2", "LAMP5", "SNCG"),
         ncol = 4,
         pt.size = -1)*theme(axis.text.x = element_text(angle = 90))
dev.off()

Idents(harmonized_object) <- harmonized_object$Sample

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/3_SexSpecific_Vlns.pdf",
         height = 8, width = 12)
DotPlot(harmonized_object, features = c("UTY", "SRY", "NLGN4Y", "XIST", "TSIX", "FIRRE"))+ theme(axis.text.x = element_text(angle = 90))
dev.off()

#prep data

Idents(harmonized_object) <- harmonized_object$RNA_snn_res.0.7 
harmonized_object <- DietSeurat(harmonized_object)
variable_genes <- VariableFeatures(harmonized_object)
harmonized_object@meta.data["study_id"] <- rep("harmonized_all", dim(harmonized_object)[2])
harmonized_object <- as.SingleCellExperiment(harmonized_object)
harmonized_object <- harmonized_object[(rownames(harmonized_object) %in% variable_genes), ]
 
#train model on variable genes
pretrained_model = MetaNeighbor::trainModel(
 var_genes = rownames(harmonized_object),
 dat = harmonized_object,
 study_id = harmonized_object$study_id,
 cell_type = harmonized_object$RNA_snn_res.0.7
)
 
#save model to file
write.table(pretrained_model, paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/", 
	"MetaNeighbor/3_pretrained_harmonized_object_RNA_snn_res.0.7.txt"))
 
 
#male data
 
rm(list = setdiff(ls(), "do_trained_MNUS")) 
 
male_data <- readRDS(file = "/home/malosree/scratch/snRNA_seq_SObj.Rds")
male_data <- DietSeurat(male_data) 
male_data@meta.data["study_id"] <-  rep("male", dim(male_data)[2])
male_data <- as.SingleCellExperiment(male_data)
 
do_trained_MNUS("male_data", "study_id", "cluster")
 
# #STAB
# #Using all datasets containing adult human cortical cells. Obtained from http://stab.comp-sysbio.org/help/
 
rm(list = setdiff(ls(), "do_trained_MNUS"))

files_to_load <- c("h5_a", "h7", "h10", "h14")
for(file in files_to_load) load(paste0("/home/malosree/scratch/STAB_data/", file, ".RData"))

list_of_objects <- list()
for(obj in ls(pattern = "^h")) list_of_objects[[obj]] <- get(obj)
rm(list = ls(pattern = "^h"))

STAB_merged <- merge(list_of_objects[[1]], list_of_objects[c(2:4)])
STAB_merged <- DietSeurat(STAB_merged)
STAB_merged@meta.data["study_id"] <-  rep("STAB", dim(STAB_merged)[2])
STAB_merged <- NormalizeData(STAB_merged, normalization.method = "LogNormalize", scale.factor = 10000)
STAB_merged <- as.SingleCellExperiment(STAB_merged, assay = "RNA")
 
do_trained_MNUS("STAB_merged", "study_id", "cluster")
#This should be cluster1 and not priCluster. Edited on 2021.12.06
do_trained_MNUS("STAB_merged", "study_id", "cluster1")

#Use Allen Brain motor cortex referece - using human M1 cortex data from https://portal.brain-map.org/atlases-and-data/rnaseq

rm(list = setdiff(ls(), "do_trained_MNUS"))

library(tidyverse)

ABI_metadata <- read_csv("/home/malosree/scratch/ABI_M1_data/metadata.csv")
ABI_metadata <- column_to_rownames(ABI_metadata, "sample_name")

ABI_data <- read_csv("/home/malosree/scratch/ABI_M1_data/matrix.csv")
ABI_data <- column_to_rownames(ABI_data, "sample_name")

ABI_metadata <- as.data.frame(ABI_metadata)
ABI_data <- t(as.matrix(ABI_data))
ABI_data <- as(ABI_data, "dgCMatrix")

ABI_data <- CreateSeuratObject(ABI_data, meta.data = ABI_metadata)
ABI_data <- NormalizeData(ABI_data, normalization.method = "LogNormalize", scale.factor = 10000)
ABI_data@meta.data["study_id"] <- rep("ABI", dim(ABI_data)[2])
ABI_data <- as.SingleCellExperiment(ABI_data)

do_trained_MNUS("ABI_data", "study_id", "cluster_label")
do_trained_MNUS("ABI_data", "study_id", "subclass_label")

print(summary(warnings()))
sessionInfo()
