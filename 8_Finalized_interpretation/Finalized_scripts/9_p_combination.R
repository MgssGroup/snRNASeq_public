library(metaRNASeq)
library(tidyverse)
library(patchwork)

base_path <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/"

directories <- c(male_broad = "04_Male_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Male_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
                male_subtype = "Mar9_2022_updated_res/",
                female_broad =  "03_Female_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
                female_subtype = "03_Female_dataset_analysis/01_pseudobulk/02_cell_subtype/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_subtype/")

types <- list(subtype = c("male_subtype", "female_subtype"), broad = c("male_broad", "female_broad"))

pdf(file = paste0(base_path, "Finalized_outputs/9_p_values_dist.pdf"), onefile = TRUE)
for(this_type_name in names(types)) {
	this_type <- types[[this_type_name]]
	male_data <-  read_csv(file = paste0(base_path, directories[[this_type[1]]], "01_Table_raw_masterfile_cpm_Freq_DS_stats_greater3.csv")) %>% 
		select(-c(1))
	female_data <-  read_csv(file = paste0(base_path, directories[[this_type[2]]], "01_Table_raw_masterfile_cpm_Freq_DS_stats_greater3.csv")) %>%
                select(-c(1))
	cell_clusters <- intersect(unique(male_data$cluster_id), unique(female_data$cluster_id))
	all_cluster_data <- list()
	all_cluster_data_unfiltered <- list()
	plots_list <- list()
	for(this_cluster in cell_clusters) {
		male_data_subset <- male_data %>% filter(cluster_id == this_cluster, Greater3 == TRUE) %>%  
			select(gene, cluster_id, logFC, p_val, p_adj.loc)
		female_data_subset <- female_data %>% filter(cluster_id == this_cluster, Greater3 == TRUE) %>%  
			select(gene, cluster_id, logFC, p_val, p_adj.loc) 
		subset_combined <- inner_join(male_data_subset, female_data_subset, by = "gene", suffix = c(".male", ".female")) 
		subset_combined %>% select(contains("p_val")) -> ind_p_val
		ind_p_val <-list( male = ind_p_val$p_val.male, female = ind_p_val$p_val.female)
		fisher_comb_res <- fishercomb(ind_p_val, BHth = 0.05)
		all_cluster_data[[this_cluster]] <- bind_cols(subset_combined, fisher_comb_res[2:4]) %>% filter(adjpval < 0.05)
		all_cluster_data_unfiltered[[this_cluster]] <- bind_cols(subset_combined, fisher_comb_res[2:4]) 
		MASS::truehist(ind_p_val$male, xlab = paste("Male p_values", this_cluster))
		MASS::truehist(ind_p_val$female, xlab = paste("Female p_values", this_cluster))
	}
	print(plots_list)		
	all_cluster_data <- bind_rows(all_cluster_data) %>% mutate(signs = sign(logFC.male)*sign(logFC.female))
	all_cluster_data_unfiltered <- bind_rows(all_cluster_data_unfiltered) %>% mutate(signs = sign(logFC.male)*sign(logFC.female))
	write.csv(all_cluster_data, file = paste0(base_path, "Finalized_outputs/9_combination_pvals_all_", this_type_name, ".csv"))
	write.csv(all_cluster_data_unfiltered, file = paste0(base_path, "Finalized_outputs/9_combination_pvals_all_unfiltered", this_type_name, ".csv"))	
}
dev.off()

print(summary(warnings()))
sessionInfo()
