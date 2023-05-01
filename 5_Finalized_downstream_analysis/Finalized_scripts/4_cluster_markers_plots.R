library(Seurat)
library(tidyverse)
library(ggplot2)

harmonized_object <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

main_markers <- c("SNAP25", "RBFOX3", "GRIN1", "SLC17A7", "GAD1", "GAD2", "ALDH1L1", "AQP4", "GFAP", "MBP", "PLP1", "PDGFRA", "SOX10", "CX3CR1", "CLDN5")

excitatory_layer_markers <- c("SLC17A7", "SATB2", "FXYD6", "GSG1L", "RASGRF2", "CUX2", "RELN", "TLE1", "RORB", "PDE1A", "PCP4", "BCL11B", "HTR2C", "RXFP1", 
		"SYNPR", "NR4A2", "TOX", "FOXP2", "TLE4", "ETV1", "NTNG2")

inhibitory_subtype_markers <- c("GAD1", "GAD2", "SST", "PVALB", "VIP", "LAMP5", "LHX6", "ADARB2", "CCK", "NPY", "CALB1", "CALB2")

glial_markers <- c("GFAP", "GJA1", "ALDH1A1", "AQP4", "MAG", "MOG", "OLIG1", "OLIG2", "SOX10", "MYT1", "PDGFRA", "ZFPM2", "ITPR2", "TCF7L2", "MRC1", "CX3CR1", "SPI1", "VIM", "CLDN5")

Idents(harmonized_object) <- harmonized_object$Cluster

all <- Reduce(union, list(main_markers, excitatory_layer_markers, inhibitory_subtype_markers, glial_markers, "nCount_RNA", "nFeature_RNA"))

inhibitory_types <- grep("^In", levels(Idents(harmonized_object)), value = TRUE)
excitatory_types <- grep("^Ex", levels(Idents(harmonized_object)), value = TRUE)
glial_types <- Reduce(union, list(grep("^Ast", levels(Idents(harmonized_object)), value = TRUE), grep("^O", levels(Idents(harmonized_object)), value = TRUE), "End", "Mic") )

feature_plots_list <- list()

violin_plots_list <- list()

heatmaps_list <- list()

i <- 1

while(i < length(all)) {

	if(i+3 <= length(all)) {
		feature_plots_list[[i]] <- FeaturePlot(harmonized_object, features = all[i:(i+3)], ncol=2, reduction = "UMAPHarmonySeededBatchSampleChemistry")
	}
	else {
		feature_plots_list[[i]] <- FeaturePlot(harmonized_object, features = all[i:length(all)], ncol = 2, reduction = "UMAPHarmonySeededBatchSampleChemistry")
	}
	i <- i+4
}

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Celltype_marker_plots/4_FeaturePlotMarkers.pdf", 
		height = 10, width = 12, onefile = TRUE)
	DimPlot(harmonized_object, reduction = "UMAPHarmonySeededBatchSampleChemistry", label = TRUE)
	feature_plots_list
dev.off()

i <- 1

while(i < length(excitatory_layer_markers)) {

        if(i+3 <= length(excitatory_layer_markers)) {
                violin_plots_list[[paste0("Ex", i)]] <- VlnPlot(harmonized_object, features = excitatory_layer_markers[i:(i+3)], ncol=2, pt.size = -1, 
									idents = excitatory_types)
        }
        else {
                violin_plots_list[[paste0("Ex", i)]] <- VlnPlot(harmonized_object, features = excitatory_layer_markers[i:length(excitatory_layer_markers)], ncol = 2, 
								pt.size = -1, idents = excitatory_types)
        }
        i <- i+4
}

i <- 1

while(i < length(inhibitory_subtype_markers)) {

        if(i+3 <= length(inhibitory_subtype_markers)) {
                violin_plots_list[[paste0("In", i)]] <- VlnPlot(harmonized_object, features = inhibitory_subtype_markers[i:(i+3)], ncol=2, pt.size = -1,
                                                                        idents = inhibitory_types)
        }
        else {
                violin_plots_list[[paste0("In", i)]] <- VlnPlot(harmonized_object, features = inhibitory_subtype_markers[i:length(inhibitory_subtype_markers)], ncol = 2,
                                                                pt.size = -1, idents = inhibitory_types)
        }
        i <- i+4
}


i <- 1

while(i < length(glial_markers)) {

        if(i+3 <= length(glial_markers)) {
                violin_plots_list[[paste0("Gl", i)]] <- VlnPlot(harmonized_object, features = glial_markers[i:(i+3)], ncol=2, pt.size = -1,
                                                                        idents = glial_types)
        }
        else {
                violin_plots_list[[paste0("Gl", i)]] <- VlnPlot(harmonized_object, features = glial_markers[i:length(glial_markers)], ncol = 2,
                                                                pt.size = -1, idents = glial_types)
        }
        i <- i+4
}

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Celltype_marker_plots/4_VlnMarkers.pdf", onefile = TRUE)
	VlnPlot(harmonized_object, features = c("nCount_RNA", "nFeature_RNA"), pt.size = -1, ncol =1)
	violin_plots_list
dev.off()

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Celltype_marker_plots/4_DotPlotsMarkers.pdf", onefile = TRUE)
	DotPlot(harmonized_object, idents = excitatory_types, features = excitatory_layer_markers)+theme(axis.text.x = element_text(angle = 90))
	DotPlot(harmonized_object, idents = inhibitory_types, features = inhibitory_subtype_markers)+theme(axis.text.x = element_text(angle = 90))
	DotPlot(harmonized_object, idents = glial_types, features = glial_markers)+theme(axis.text.x = element_text(angle = 90))
	DotPlot(harmonized_object, features = main_markers)+theme(axis.text.x = element_text(angle = 90))
dev.off()

#number of UMIs and genes and cells per cluster

harmonized_object@meta.data %>% group_by(Cluster) %>% select(nCount_RNA, nFeature_RNA) %>% 
		summarise(mean_UMI = mean(nCount_RNA),
				median_UMI = median(nCount_RNA),
				mean_genes = mean(nFeature_RNA), 
				median_genes = median(nFeature_RNA),
				cells = n()) -> count_summary_all

harmonized_object@meta.data %>% group_by(Cluster, Chemistry) %>% select(nCount_RNA, nFeature_RNA) %>%
                summarise(mean_UMI = mean(nCount_RNA),
                                median_UMI = median(nCount_RNA),
                                mean_genes = mean(nFeature_RNA),
                                median_genes = median(nFeature_RNA),
                                cells = n()) -> count_summary_chemistry

harmonized_object@meta.data %>% group_by(Cluster, Sex) %>% select(nCount_RNA, nFeature_RNA) %>%
                summarise(mean_UMI = mean(nCount_RNA),
                                median_UMI = median(nCount_RNA),
                                mean_genes = mean(nFeature_RNA),
                                median_genes = median(nFeature_RNA),
                                cells = n()) -> count_summary_dataset

write.csv(count_summary_all, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Celltype_marker_plots/4_count_summary_all.csv")
write.csv(count_summary_chemistry, 
	file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Celltype_marker_plots/4_count_summary_chemistry.csv")
write.csv(count_summary_dataset, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Celltype_marker_plots/4_count_summary_dataset.csv")

print(summary(warnings()))
sessionInfo()
