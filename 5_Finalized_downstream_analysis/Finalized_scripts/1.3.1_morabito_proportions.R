#Inspired by https://github.com/swaruplab/Single-nuclei-epigenomic-and-transcriptomic-landscape-in-Alzheimer-disease/blob/master/snATAC_cluster_characterization.Rmd
#Morabito et al., 2021 (Swarup lab)

library(Seurat)
library(tidyverse)
library(ggplot2)
library(rstatix)

#set.seed(21)

#Load data
print("Load data")
harmonized_all <- readRDS("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

harmonized_all$Sample <- str_replace(harmonized_all$Sample, pattern = fixed("M24_2"), replacement = "M24")

cluster_cats <- c("Cluster", "Broad")
sexes <- list(c("Male", "Female"), c("Male"), c("Female"))
sample_percent <- 0.70
n_iterations  <- 100

for(sex in sexes) {

set.seed(21)

harmonized_object <- subset(harmonized_all, Sex %in% sex)

for (cluster_cat in cluster_cats) {

cur_sample <- harmonized_object@meta.data

#Redo baseline calculations

#Get subject proportions in each cluster
print("Get subject proportions for each cluster using Morabito method")
meta_list <- cur_sample %>% dplyr::group_split(Sample)
temp <- lapply(meta_list, function(meta){
	        df <- as.data.frame(meta[,cluster_cat] %>% table / nrow(meta))
		colnames(df) <- c('cluster', 'proportion')
                df$Sample <- unique(meta$Sample)
                df$Condition <- unique(meta$Condition)
                df
        })

cur_df <- Reduce(rbind, temp)

clusters <- unique(as.character(t(harmonized_object@meta.data[cluster_cat])))

print("Recalculate baseline results to make sure they are the same as those with my method.")
#This does not seem to give the same numbers when it is re-run multiple times
lapply(clusters, function(cur_cluster){
    print(cur_cluster)
    result <- wilcox_test(
	cur_df[cur_df$cluster == cur_cluster,],
	proportion~Condition,
	paired = FALSE,
	alternative = "two.sided",
	detailed = TRUE	
  )
  result$cluster  <- cur_cluster
  result
}) %>% bind_rows() %>% mutate(p_adjust = p.adjust(p, method = "BH"))-> baseline_results

sex <- paste0(sex, collapse ="")

write.csv(baseline_results, file = paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Morabito_scDC/1.3.1_baseline_results_",
				cluster_cat, "_", sex,
				".csv"))

proportion_df <- data.frame()
p_values <- data.frame()

for (i in 1:n_iterations) {

	#Subsample Seurat object 
	print("Subsample Seurat object metadata")
	cur_sample <- harmonized_object@meta.data[sample(rownames(harmonized_object@meta.data), round(sample_percent*ncol(harmonized_object))),]

	#Get subject proportions in each cluster
	print("Get subject proportions for each cluster using subsampled metadata")
	meta_list <- cur_sample %>% dplyr::group_split(Sample)
	temp <- lapply(meta_list, function(meta){
		df <- as.data.frame(meta[,cluster_cat] %>% table / nrow(meta))
		colnames(df) <- c('cluster', 'proportion')
		df$Sample <- unique(meta$Sample)
		df$Condition <- unique(meta$Condition)
		df
	})

	cur_df <- Reduce(rbind, temp)
	cur_df$iteration <- i

	clusters <- unique(as.character(t(cur_sample[cluster_cat])))
	
	#This is a custom addition
	print(c("Do the Wilxocon test from the proportions for the subsetted data for each cluster", i))
	lapply(clusters, function(cur_cluster){
        #print(cur_cluster)
        result <- wilcox_test(
        cur_df[cur_df$cluster == cur_cluster,],
        proportion~Condition,
        paired = FALSE,
        alternative = "two.sided",
        detailed = TRUE)
	result$cluster <- cur_cluster 
	result$i <- i
        result
        }) %>% bind_rows() -> wilcox_res
	
	#This is a custom addition
	print("Store the Wilcoxon results from this iteration of subsetting")
	p_values <- rbind(p_values, wilcox_res)	

	proportion_df <- rbind(proportion_df, cur_df)

}

#This is a custom addition	
p_values %>% group_by(cluster) %>% 
	group_map(~list(Cluster = unique(.$cluster), 
			booted_p_min = min(.$p), 
			booted_p_med = median(.$p), 
			booted_p_max = max(.$p)), .keep = TRUE) %>% 
		bind_rows() -> booted_p_values

write.csv(booted_p_values, file = paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/", 
	"Morabito_scDC/1.3.1_subsampled_booted_p", 
	cluster_cat, "_", sex,
	".csv"))

clusters <- unique(as.character(t(harmonized_object@meta.data[cluster_cat])))

print("Redo Wilcoxon tests using compiled proportions from subsampling of cells to get bootstrapped p values.")
lapply(clusters, function(cur_cluster){
	print(cur_cluster)
	result <- wilcox_test(
        proportion_df[proportion_df$cluster == cur_cluster,],
	proportion~Condition,
	paired = FALSE,
	alternative = "two.sided",
	detailed = TRUE)  
	result$cluster <- cur_cluster
	result
	}) %>% bind_rows() %>% mutate(p_adjust = p.adjust(p, method = "BH")) -> Morabito_bootstrapped_results

write.csv(Morabito_bootstrapped_results, file = paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Morabito_scDC/1.3.1_",
                cluster_cat, "_", sex,
                "_Morabito_bootstrapped_results.csv"))

}
}

print(summary(warnings()))
sessionInfo()
