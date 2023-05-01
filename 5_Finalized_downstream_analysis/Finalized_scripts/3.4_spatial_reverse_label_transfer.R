#Reverse label transfer from Lieber spatial data

library(Matrix)
library(Seurat)
library(SeuratData)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(scclusteval)
library(pheatmap)
library(grid)

#Load the single-cell object and SCTransform
print("Loading the single cell object and marker information.")
harmonized_object <- readRDS("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

print("Running SCTransrform and PCA.")
harmonized_object <- SCTransform(harmonized_object, vars.to.regress = c("percent.mt", "nCount_RNA", "Chemistry", "Batch", "Sample"),
                        conserve.memory = TRUE, verbose = TRUE) %>%
                        RunPCA(npcs = 100, verbose = TRUE)

sections <- c("151673", "151507")

plots_list <- list()

for(section in sections) {

        print("Get the layer information from the spatialLIBD package.")
        load(paste0("/home/malosree/scratch/Maynard_data/section_", section, ".Rda"))
        rm(list = paste0(c("coords_", "counts_"), section))

        print("Get the 10X output from the AWS download.")
        data_dir <- paste0("/home/malosree/scratch/Maynard_data/", section)
        filename <- paste0(section, "_filtered_feature_bc_matrix.h5")
        section_Seu <- Load10X_Spatial(data_dir, filename)
        section_Seu <- AddMetaData(section_Seu, get(paste0("layers_", section)), col.name = "Layer_short")

        print("Running SCTransrform and PCA.")
        section_Seu <- SCTransform(section_Seu, assay = "Spatial", verbose = TRUE, conserve.memory = TRUE) %>%
                        RunPCA(npcs = 50, versbose = TRUE)

        print("Find transfer anchors.")
        anchors <- FindTransferAnchors(query = harmonized_object, reference = section_Seu, normalization.method = "SCT", reduction = "cca")
	
	print("Transfer data (predict layer labels for snRNA-seq clusters)")
        Layer_predictions <- TransferData(anchorset = anchors, refdata = section_Seu$Layer_short, prediction.assay = FALSE,
		weight.reduction = harmonized_object[["pca"]], dims = 1:50)
	print("Adding prediction info to MetaData")
        harmonized_object <- AddMetaData(harmonized_object, Layer_predictions)
	
	write.csv(Layer_predictions, file = paste0(
                "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/3.4_",
                section,
                "_Layer__preds.csv"))
        plots_list[[paste0(section,"_section_UMAP")]] <- DimPlot(harmonized_object, group.by = "predicted.id", 
					reduction = "UMAPHarmonySeededBatchSampleChemistry")

}

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/3.4_layer_prediction_UMAP.pdf",
	onefile = TRUE)
	DimPlot(harmonized_object, group.by = "Cluster",
                                        reduction = "UMAPHarmonySeededBatchSampleChemistry", label = TRUE) + theme(legend.position = "none")
	plots_list
	predictions <- read_csv(file = paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/3.4_",
                "151507",
                "_Layer__preds.csv"))
	harmonized_object$predicted.id <- predictions$predicted.id
	JaccardSets_matrix <- PairWiseJaccardSets(harmonized_object$Cluster, harmonized_object$predicted.id)
        pheatmap(JaccardSets_matrix, cluster_rows = FALSE, cluster_cols = FALSE)
	grid.newpage()
	predictions <- read_csv(file = paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/3.4_",
                "151673",
                "_Layer__preds.csv"))
        harmonized_object$predicted.id <- predictions$predicted.id
        JaccardSets_matrix <- PairWiseJaccardSets(harmonized_object$Cluster, harmonized_object$predicted.id)
        pheatmap(JaccardSets_matrix, cluster_rows = FALSE, cluster_cols = FALSE)

dev.off()

print(summary(warnings()))
sessionInfo()
