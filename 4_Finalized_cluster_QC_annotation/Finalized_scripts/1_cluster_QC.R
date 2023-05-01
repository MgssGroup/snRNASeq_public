library(Seurat)
library(readr)
library(dplyr)
library(ggplot2)

#########
#Pure QC#
#########

#Stability comparison to the scclusteval clusters

print("Loading scclusteval data")
subsample_idents_list <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/subsample_idents_list.Rds")

subsample_idents_list %>% filter(resolution == 0.7, k_param == 30, pc == 70) %>% 
		select(stable_cluster_80_0.8, stable_cluster_70_0.7, stable_cluster_boot) -> subsample_idents_list

fullsample_idents <- readRDS(file =
  "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/gather_full_sample.rds")

fullsample_idents %>%
		mutate(cluster_num = purrr::map_dbl(original_ident_full, ~n_distinct(.x))) %>% 
		filter(resolution == 0.7, k_param == 30, pc == 70) -> fullsample_idents

#Add the scclusteval evaluated cluster labels 
scclusteval_cluster_labels <- tibble(fullsample_idents$original_ident_full[[1]])
colnames(scclusteval_cluster_labels) <- "scclusteval_0.7_30_70"
rownames(scclusteval_cluster_labels) <- names(fullsample_idents$original_ident_full)

print("Loading Seurat object")
harmonized_object <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/5_reharmonized_seeded.Rds")

#using most lenient definition
scclusteval_cluster_labels %>% mutate(stable_70_0.7_scclusteval = plyr::revalue(scclusteval_0.7_30_70, 
	subsample_idents_list$stable_cluster_70_0.7[[1]]$stable_cluster)) -> scclusteval_cluster_labels

harmonized_object <- AddMetaData(harmonized_object, scclusteval_cluster_labels)

#find cluster and stability equivalencies for the clustering after seeding
JaccardSets_matrix <- scclusteval::PairWiseJaccardSets(unlist(fullsample_idents$original_ident_full), Idents(harmonized_object))

#just a check
dim(JaccardSets_matrix)

#cluster equivalencies
cluster_equivalencies <-  tibble(seeded_cluster = 0:40, sccluseval_clusters = apply(JaccardSets_matrix, 2, function(x) {
         which(x == max(x), arr.ind = TRUE)-1}))
cluster_equivalencies %>% mutate(stability_equivalence_70_0.7 = plyr::revalue(
	as.character(sccluseval_clusters), subsample_idents_list$stable_cluster_70_0.7[[1]]$stable_cluster)) -> cluster_equivalencies
harmonized_object@meta.data$stable_equivalent_70_0.7 <- plyr::mapvalues(Idents(harmonized_object), cluster_equivalencies$seeded_cluster, 
	cluster_equivalencies$stability_equivalence_70_0.7)


#Doublet Assessment

print("Load doublet data")
 
doublets_predictions <- read.csv(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/Finalized_outputs/6_doublet_predictions.csv", 
	row.names =1)

harmonized_object <- AddMetaData(harmonized_object, metadata = doublets_predictions)

rm(doublets_predictions)

#plot some UMAPS

print("Plot UMAPs")

plots_list <- list()

variables <- c("RNA_snn_res.0.7", "DF.classifications", "stable_equivalent_70_0.7")
for (variable in variables) {
	Idents(harmonized_object) <- harmonized_object@meta.data[,variable]
        plots_list[[variable]] <- DimPlot(harmonized_object, reduction = "UMAPHarmonySeededBatchSampleChemistry", label = TRUE)
}

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/1_cluster_QC_UMAP_plots.pdf", onefile = TRUE, height = 10, width = 12)
plots_list
dev.off()

#plot some Vlns

#Plot some basic QC parameters

print("Plot QC violins")
plots_list <- list()

features_to_plot <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "pANN")
idents_to_use <- list(all = c("RNA_snn_res.0.7"), by_sex = c("RNA_snn_res.0.7", "Sex"))
log_scale <- c(TRUE, FALSE)

for(idents in idents_to_use) {
	if (length(idents) == 1) Idents(harmonized_object) <- harmonized_object@meta.data[, idents]
	else Idents(harmonized_object) <- paste(harmonized_object@meta.data[,idents[1]], harmonized_object@meta.data[,idents[2]])
	for(scale in log_scale) {
		plots_list[[paste(c(scale, idents), collapse  = "_")]] <- VlnPlot(harmonized_object, features = features_to_plot,
		log = scale,
		ncol = 1,
		pt.size=-1)*theme(axis.text.x = element_text(angle = 90))
	}
}

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/1_per_cluster_nCount_nFeature_pt.mito_pt.ribo_DF.pdf",
        onefile = TRUE, height = 13, width = 13)
plots_list
dev.off()

#Subject, cluster, batch, chemistry, condition contributions
print("Prep to plot heatmaps")

rm(list = setdiff(ls(), "harmonized_object"))

plots_list <- list()

variables <- c("Sample", "Batch", "Chemistry", "Condition", "Sex", "DF.classifications")

scale_dirs <- c("row", "column")

subset_to_use <- list(male = c("Male"), female = c("Female"), all = c("Male", "Female")) 

Idents(harmonized_object) <- harmonized_object@meta.data[,"RNA_snn_res.0.7"]

print("Plotting heatmaps")
for(subset_var_name in names(subset_to_use)) {
	if("plotting_subset" %in% ls()) rm(list = c("plotting_subset"))
	print("subsetting")
	subset_var <- subset_to_use[[subset_var_name]]
	plotting_subset <- subset(harmonized_object, subset = Sex %in% subset_var)
	for(variable in variables) {
		print("counting")
		frequencies <- table(plotting_subset@meta.data[,variable], Idents(plotting_subset))
                write.csv(frequencies, file = paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/", 
				subset_var_name,"/1_", variable, "_", subset_var_name, "_cluster_freq.csv"))
                perc_rep <- colSums(frequencies>0)/dim(frequencies)[1]
                write.csv(perc_rep, file = paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/", 
				subset_var_name,"/1_", variable, "_", subset_var_name, "_cluster_perc_rep.csv"))
		print("plotting")
		if(!(variable %in% c("Chemistry", "Sex") & subset_var_name == "male") & !(variable == "Sex" & subset_var_name == "female")) {
			pdf(file = paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/", 
				subset_var_name, "/1_",variable, "_cluster_categorical_counts_", subset_var_name,".pdf"),
				onefile = TRUE, height = 13, width = 13)
			for(scale_dir in scale_dirs) {
				pheatmap::pheatmap(frequencies[,colnames(frequencies)[order(as.numeric(colnames(frequencies)))]], 
					cluster_cols = FALSE, cluster_rows = FALSE, scale = scale_dir)	
			}
			dev.off()
		}
	}
}  

#PCA of cell-type proportions for each subject by coloured by condition, batch, chemistry, sex

print("Subject cell-type proportion PCA")

rm(list = setdiff(ls(), "harmonized_object"))

subset_to_use <- list(male = c("Male"), female = c("Female"), all = c("Male", "Female"))

plots_list <- list()

for(subset_var_name in names(subset_to_use)) {
	if("plotting_subset" %in% ls()) rm(list = c("plotting_subset"))
        print("subsetting")
        subset_var <- subset_to_use[[subset_var_name]]	
	plotting_subset <- subset(harmonized_object, subset = Sex %in% subset_var)
	subject_cluster_proportions <- table(plotting_subset$Sample, Idents(plotting_subset))/rowSums(table(plotting_subset$Sample, 
					Idents(plotting_subset)))

	print("Get Batch, Chemistry, Sex, Condition info in the correct order")
	#The order was checked in this case and matched because the metadata subject names were ascending and alphabetical, but this would not work in the general 
	#case as the ouput of table orders things alphabetically.

	data_for_pca <- data.frame(Subject = unique(plotting_subset$Sample),
		Condition = stringr:: str_split(unique(paste(plotting_subset$Sample, plotting_subset$Condition)), pattern = " ", simplify = TRUE)[,2],
		Batch =  stringr:: str_split(unique(paste(plotting_subset$Sample, plotting_subset$Batch)), pattern = " ", simplify = TRUE)[,2],
		Chemistry =  stringr:: str_split(unique(paste(plotting_subset$Sample, plotting_subset$Chemistry)), pattern = " ", simplify = TRUE)[,2],
		Sex =  stringr:: str_split(unique(paste(plotting_subset$Sample, plotting_subset$Sex)), pattern = " ", simplify = TRUE)[,2])
	
	#check if direction of this pca makes sense 
	print("Running PCA and making plots") 
	res_pca <- prcomp(subject_cluster_proportions, scale. = TRUE)
	colour_vars <- c("Condition", "Batch", "Chemistry")
	for(colour_var in colour_vars) {
		plots_list[[paste0(colour_var, subset_var_name)]] <- ggbiplot::ggbiplot(res_pca, groups = data_for_pca[,colour_var], ellipse = TRUE, var.axes = FALSE)
	}
	if(subset_var_name == "all") {
		plots_list[[paste0("Sex", subset_var_name)]] <- ggbiplot::ggbiplot(res_pca, groups = data_for_pca$Sex, ellipse = TRUE, var.axes = FALSE) 
	}
} 

pdf(file = paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/1_Subject_cluster_proportions_PCA_.pdf"),
                   onefile = TRUE, height = 10, width = 12)
plots_list
dev.off()

#PCA of subject proportions for each cell-type

print("Cell-type subject proportion PCA")

rm(list = setdiff(ls(), "harmonized_object"))

subset_to_use <- list(male = c("Male"), female = c("Female"), all = c("Male", "Female"))

plots_list <- list()

for(subset_var_name in names(subset_to_use)) {
        if("plotting_subset" %in% ls()) rm(list = c("plotting_subset"))
        print("subsetting")
        subset_var <- subset_to_use[[subset_var_name]]
        plotting_subset <- subset(harmonized_object, subset = Sex %in% subset_var)
        cluster_subject_proportions <- table(plotting_subset$RNA_snn_res.0.7, plotting_subset$Sample)/rowSums(table(plotting_subset$RNA_snn_res.0.7, plotting_subset$Sample))
        #check if direction of this pca makes sense
        print("Running PCA and making plots")
        res_pca <- prcomp(na.omit(cluster_subject_proportions), scale. = TRUE)
	data_for_pca <- cbind(as.data.frame(res_pca$x), clusters = rownames(res_pca$x))
        plots_list[[subset_var_name]] <- ggplot(data = data_for_pca, aes(x=PC1, y = PC2, label = clusters)) + geom_text() 
}

pdf(file = paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/1_Cluster_subject_proportions_PCA_.pdf"),
                   onefile = TRUE, height = 10, width = 12)
plots_list
dev.off()


#Repeat the above PCA analysis again with the retained clusters

saveRDS(harmonized_object, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

print(summary(warnings()))
sessionInfo()

