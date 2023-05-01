library(tidyverse)
library(patchwork)
library(muscat)
library(scater)

base_path <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/"

directories <- c("male_broad"=  "04_Male_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Male_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
                "male_subtype" = "Mar9_2022_updated_res/",
                "female_broad" = "03_Female_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
                "female_subtype" = "03_Female_dataset_analysis/01_pseudobulk/02_cell_subtype/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_subtype/")

metadata <- read_csv(file = paste0(base_path, "Male_female_metadata_combined.csv")) %>% mutate(Sample = paste0(Sample, ".cpm"))

print("Calculate the number of cells per cluster and broad cell-type for each sex (this will included the cells from subjects with low contribution)")
harmonized_object <- readRDS(file = paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds"))

Sexes <- list("male" = "Male", "female" = "Female")
Categories <- list("broad" = "Broad", "subtype" = "Cluster")

cell_numbers <- list()

for(sex in names(Sexes)) {
	for(category in names(Categories)) {
		res <- table(
			harmonized_object@meta.data[[Categories[[category]]]][harmonized_object@meta.data[["Sex"]]==Sexes[[sex]]])
		res <- enframe(res, name = "cluster", value = "num_cells")
		cell_numbers[[paste0(sex, "_", category)]] <- res	
	}
}

find_condition <- function(subject_name) {
        condition <- metadata$Condition[metadata$Sample == subject_name]
        condition
}

print("Make per cluster boxplots showing the CPM for each subject for each DEG.")

for(directory_name in names(directories)) {
	print(directory_name)
	directory <- directories[[directory_name]]
	
	filtered_results <- read_csv(file = paste0(base_path, directory, "01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct.csv")) %>%
                     mutate(cluster_id = recode(cluster_id, InN8_ADARB2 = "InN8_Mix"))
	
	pdf(file = paste0(base_path, "Finalized_outputs/Diff_exp_plots/Subject_Boxplots/1.7_", directory_name, "_distribution_zeros.pdf"))
		par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
		p1 <- MASS::truehist(filtered_results$NumNonZero, h = 1, prob = FALSE, xlab = "NumNonZero")
		print(p1)
		write.csv(filtered_results$NumNonZero, file = 
			paste0(base_path, "Finalized_outputs/Diff_exp_plots/Subject_Boxplots/1.7_source_data_",  directory_name, "_distribution_zeros.csv"))	
		print(MASS::truehist(filtered_results$NumNonZero[filtered_results$cluster_id != "Mix"], h =1, prob = FALSE), xlab = "NumNonZero")
	dev.off()
	
	#Log version added 2022.04.18
	pdf(file = paste0(base_path, "Finalized_outputs/Diff_exp_plots/Subject_Boxplots/1.7_", directory_name, "_distribution_outliers.pdf"))
		par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
                p1 <- MASS::truehist(filtered_results$NumOutliers, h = 1, prob = FALSE, xlab = "NumOutliers")
		print(p1)
                write.csv(filtered_results$NumOutliers, file =
                        paste0(base_path, "Finalized_outputs/Diff_exp_plots/Subject_Boxplots/1.7_source_data_",  directory_name, "_distribution_outliers.csv"))
                print(MASS::truehist(filtered_results$NumOutliers[filtered_results$cluster_id != "Mix"], h = 1, prob = FALSE, xlab = "NumOutliers"))
		print(MASS::truehist(filtered_results$NumOutliersLog, h = 1, prob = FALSE, xlab = "NumOutliersLog"))
                print(MASS::truehist(filtered_results$NumOutliersLog[filtered_results$cluster_id != "Mix"], h = 1, prob = FALSE, xlab = "NumOutliersLog"))
        dev.off()
	
	pdf(file = paste0(base_path, "Finalized_outputs/Diff_exp_plots/Subject_Boxplots/1.7_", directory_name, "_num_DEG_num_sub.pdf"))
		filtered_results %>% mutate(NumSub = rowSums(across(ends_with(".cpm"), ~(.x>=0)), na.rm = TRUE)) %>% group_by(cluster_id) %>%
			summarize(num_DEG = n(), num_sub = unique(NumSub), cluster = unique(cluster_id)) -> to_plot
		p1 <- ggplot(to_plot, aes(x = num_sub, y = num_DEG, label = cluster)) + geom_label(size = 7) + theme_classic() + ggtitle(directory_name) +
			theme_classic(base_size = 22) + theme(axis.text.x = element_text(size = 22, angle = 45, hjust = 1))
		print(p1)
                write.csv(p1$data, file =
                        paste0(base_path, "Finalized_outputs/Diff_exp_plots/Subject_Boxplots/1.7_source_data_",  directory_name, "_num_DEGs_num_sub.csv"))
	dev.off()

	pdf(file = paste0(base_path, "Finalized_outputs/Diff_exp_plots/Subject_Boxplots/1.7_", directory_name, "_num_DEG_num_cells.pdf"))
                filtered_results %>% group_by(cluster_id) %>%
                        summarize(num_DEG = n(), cluster = unique(cluster_id)) -> to_plot
		to_plot <- inner_join(cell_numbers[[directory_name]], to_plot, by = "cluster")
                p1 <- ggplot(to_plot, aes(x = num_cells, y = num_DEG, label = cluster)) + geom_label(size = 7) + theme_classic() + ggtitle(directory_name) +
			theme_classic(base_size = 22) + theme(axis.text.x = element_text(size = 22, angle = 45, hjust = 1))
		print(p1)
                write.csv(p1$data, file =
                        paste0(base_path, "Finalized_outputs/Diff_exp_plots/Subject_Boxplots/1.7_source_data_",  directory_name, "_num_DEGs_num_cells.csv"))
        dev.off()

	lapply(unique(filtered_results$cluster_id), function(cluster) {
		print(cluster)
	
		#Get subjects which are not NA for this cluster
                filtered_results %>%  filter(cluster_id == cluster) %>%
                        select(any_of(c(metadata$Sample, "gene", "logFC"))) %>%
                        select(where(~!all(is.na(.x)))) %>% pivot_longer(any_of(c(metadata$Sample)), names_to = "Subject", values_to = "CPM") %>%
                        mutate(Condition = map_chr(Subject, find_condition)) %>% arrange(desc(logFC)) -> to_plot

		#Get siginificant genes to plot
                filtered_results %>%  filter(cluster_id == cluster) %>% arrange(desc(logFC)) %>% select(gene) %>% unlist() -> genes_to_plot
                genes_to_plot %>% split(ceiling(seq_along(genes_to_plot)/6)) -> genes_to_plot

                pdf(file = paste0(base_path, "Finalized_outputs/Diff_exp_plots/Subject_Boxplots/1.7_", directory_name, "_", cluster, "_boxplots.pdf"), 
			onefile = TRUE,
                        width = 18, height = 8)
                for(genes in genes_to_plot) {
			list_of_plots <- list()
			for(this_gene in genes) {
				to_plot_now <- filter(to_plot, gene == this_gene) 
				list_of_plots[[this_gene]] <- ggplot(to_plot_now, aes(x = Condition, y = CPM, label = Subject, fill = Condition)) + 
					geom_boxplot() + geom_label() + theme_classic() + 
					ggtitle(this_gene)
			}
			patched <- Reduce("+", list_of_plots) + plot_layout(ncol = 3, byrow = TRUE)
			print(patched)						
                }
                dev.off()
                NULL
        })
}

print(summary(warnings()))
sessionInfo()
