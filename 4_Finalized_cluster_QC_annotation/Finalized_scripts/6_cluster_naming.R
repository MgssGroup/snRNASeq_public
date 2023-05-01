#Step2 had to be redone right before this because of a problem with saving the all_markers file, which wiped it out 
#Not re-run again on 2021.12.07 although Step3 had to be redone to fix an error with STAB cell type labels used. 
library(Seurat)
library(tidyverse)

print("Loading data")
harmonized_object <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")
clusters_names <- read_csv(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Cluster_correspondence.csv")

clusters_names %>% mutate(original = as.character(RNA_snn_res.0.7)) %>% select(original, Cluster) %>% deframe() -> fine_clusters
clusters_names %>% mutate(original = as.character(RNA_snn_res.0.7)) %>% select(original, Broad) %>% deframe() -> broad_clusters

harmonized_object@meta.data["Cluster"] <- recode(harmonized_object$RNA_snn_res.0.7, !!!fine_clusters)
harmonized_object@meta.data["Broad"] <- recode(harmonized_object$RNA_snn_res.0.7, !!!broad_clusters)

Idents(harmonized_object) <- harmonized_object$Cluster

print("Rewriting the cluster markers file produced previously in this series of analysis")
all_markers <- read_csv(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/2_all_markers_Presto_WilcoxAUC.csv")
all_markers %>% mutate(group_named = recode(group, !!!fine_clusters)) -> all_markers

write.csv(all_markers, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/6_all_markers_Presto_WilcoxAUC.csv")

harmonized_object <- BuildClusterTree(harmonized_object, reduction = "Harmony_Seeded_Batch_Sample_Chemistry")

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/6_labelled_plots.pdf", onefile = TRUE, 
	height = 12, width = 12)
DotPlot(harmonized_object, features = c("SNAP25", "SLC17A7", "GAD1", "ALDH1L1", "PDGFRA", "PLP1", "CLDN5", "CX3CR1"))
PlotClusterTree(harmonized_object)
dev.off()		

print("Rewriting the Seurat object corresponding to this series of analysis")
saveRDS(harmonized_object, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

print(summary(warnings()))
sessionInfo()
