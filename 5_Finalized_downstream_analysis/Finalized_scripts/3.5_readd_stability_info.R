#Re-add missed stability info

library(Seurat)
library(readr)
library(dplyr)
library(ggplot2)
library(scclusteval)
library(pheatmap)

#Re-run stability comparison to the scclusteval clusters

print("Loading scclusteval data")
subsample_idents_list <- readRDS(file = paste0("/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/",
				"subsample_idents_list.Rds"))

subsample_idents_list %>% filter(resolution == 0.7, k_param == 30, pc == 70) %>%
                select(data, stable_cluster_80_0.8, stable_cluster_70_0.7, stable_cluster_boot) -> subsample_idents_list

fullsample_idents <- readRDS(file =
  "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/gather_full_sample.rds")

fullsample_idents %>%
                mutate(cluster_num = purrr::map_dbl(original_ident_full, ~n_distinct(.x))) %>%
                filter(resolution == 0.7, k_param == 30, pc == 70) -> fullsample_idents

#Add the scclusteval evaluated cluster labels
scclusteval_cluster_labels <- tibble(fullsample_idents$original_ident_full[[1]])
colnames(scclusteval_cluster_labels) <- "scclusteval_0.7_30_70"

print("Loading Seurat object")
harmonized_object <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

#using most lenient definition
scclusteval_cluster_labels %>% mutate(stable_70_0.7_scclusteval = plyr::revalue(scclusteval_0.7_30_70,
        subsample_idents_list$stable_cluster_70_0.7[[1]]$stable_cluster)) -> scclusteval_cluster_labels
rownames(scclusteval_cluster_labels) <- names(fullsample_idents$original_ident_full[[1]])

harmonized_object@meta.data$scclusteval_0.7_30_70 <- NULL
harmonized_object@meta.data$stable_70_0.7_scclusteval <- NULL

harmonized_object <- AddMetaData(harmonized_object, scclusteval_cluster_labels)

#JaccardRainCloudPlot produced from the scclusteval_round2/pyflow_seurat_parameter_custom/post_snakemake.R was for the wrong set of parameters 
#(calculations of bootstrapped median should still be correct)
pdf(file = "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/BestJaccard_correct.pdf",
        height = 10, width = 12, onefile = TRUE)
JaccardRainCloudPlot(subsample_idents_list$data[[1]]$original_ident,
                          subsample_idents_list$data[[1]]$recluster_ident) +
        geom_hline(yintercept = c(0.7, 0.8), linetype = 2) +
        xlab("cluster id w/ k=30 res=0.7 pc=70")
JaccardSets_matrix <- PairWiseJaccardSets(unlist(fullsample_idents$original_ident_full), harmonized_object$Cluster)
pheatmap(JaccardSets_matrix, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

harmonized_object <- saveRDS(harmonized_object, 
		file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

print(summary(warnings()))
sessionInfo()
