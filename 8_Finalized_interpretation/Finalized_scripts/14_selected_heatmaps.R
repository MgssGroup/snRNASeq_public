library(muscat)
library(ComplexHeatmap)
library(tidyverse)

base_path <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/"

directories <- c("male_broad"=  "04_Male_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Male_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
		"male_subtype" = "Mar9_2022_updated_res/",
		"female_broad" = "03_Female_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
		"female_subtype" = "03_Female_dataset_analysis/01_pseudobulk/02_cell_subtype/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_subtype/")

metadata <- read_csv(file = paste0(base_path, "Male_female_metadata_combined.csv")) %>% mutate(Sample = paste0(Sample, ".cpm"))

case_subjects <- metadata %>% filter(Condition == "Case") %>% select(Sample) %>% unlist()
print(case_subjects)
control_subjects <- metadata %>% filter(Condition == "Control") %>% select(Sample) %>% unlist()
print(control_subjects)

directory <- directories[["female_subtype"]]

filtered_results <- read_csv(file = paste0(base_path, directory, "/01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct.csv"))

heatmap_maker <- function(exprs_vals, cluster_id, marker = FALSE, h_dim = 0.3, f_size = 10, selected_labels = 0) {
#Remove NA subjects and convert the expression to matrix
exprs_vals %>% column_to_rownames("gene") %>% select(where(~!all(is.na(.x)))) %>% as.matrix() -> exprs_vals
		if(marker == FALSE) {
		#Scale each row (expression of each gene)
                        exprs_dim_names <- dimnames(exprs_vals)
                        exprs_vals <- t(apply(exprs_vals, 1, scale))
                        dimnames(exprs_vals) <- exprs_dim_names
                } else {
                        #Scale each col (expression of marker genes within a subject)
                        exprs_dim_names <- dimnames(exprs_vals)
                        exprs_vals <- apply(exprs_vals, 2, scale)
                        dimnames(exprs_vals) <- exprs_dim_names
                }
                #Set colors for plotting
                col_fun = circlize::colorRamp2(c(min(exprs_vals), 0, max(exprs_vals)),
                        c('#edf8b1','#7fcdbb','#2c7fb8'))
                #Create annotation data
                column_annot <- metadata %>% select(Sample, Condition) %>% filter(Sample %in% colnames(exprs_vals)) %>%
			mutate(Sample = str_sub(Sample, end = -5)) %>%
                        column_to_rownames("Sample") %>% as.data.frame()
		colnames(exprs_vals) <- str_sub(colnames(exprs_vals), end = -5)
                column_annot <- HeatmapAnnotation(df = column_annot, col = list(Condition = c("Case" = '#f1a340', "Control" = '#998ec3')))
                #Plot heatmap
                print(dim(exprs_vals)[1])
                print(dim(exprs_vals)[2])
                print(cluster_id)
		
		hm <- Heatmap(exprs_vals, column_title = cluster_id, top_annotation = column_annot, cluster_columns = FALSE,
                                row_order = seq_len(dim(exprs_vals)[1]),
                                col = col_fun,
                                height = unit(h_dim*dim(exprs_vals)[1], "cm"),
                                width = unit(0.3*dim(exprs_vals)[2], "cm"),
                                row_names_gp = grid::gpar(fontsize = f_size),
                                column_names_gp = grid::gpar(fontsize = 10))

		if(selected_labels > 0) {
			
			write.csv(rownames(exprs_vals), file = paste0(base_path, "Finalized_outputs/14_ExN10_gene_order.csv"))
			at <- c(1:selected_labels, (dim(exprs_vals)[1]-(selected_labels-1)):dim(exprs_vals)[1])
			labels <- rownames(exprs_vals)[at]			
			ra <- rowAnnotation(label = anno_mark(at = at, labels = labels))
			hm <- hm+ra

                }

		hm 

}

mic_data <- filtered_results %>% filter(cluster_id == "Mic1") %>% arrange(logFC) %>%
                                select(any_of(c(metadata$Sample, "gene")))
write.csv(mic_data, file = paste0(base_path, "Finalized_outputs/14_mic_heatmap_source_data.csv"))
svg(file = paste0(base_path, "Finalized_outputs/14_Female_Mic1_heatmap.svg"), height = 10, width = 8)
	heatmap_maker(mic_data, "Mic1", marker = FALSE) %>% draw
dev.off()

InN1_PV_data <- filtered_results %>% filter(cluster_id == "InN1_PV") %>% arrange(logFC) %>%
                                select(any_of(c(metadata$Sample, "gene")))
write.csv(InN1_PV_data, file = paste0(base_path, "Finalized_outputs/14_InN1_PV_heatmap_source_data.csv"))
svg(file = paste0(base_path, "Finalized_outputs/14_Female_InN1_PV_heatmap.svg"), height = 7, width = 9)
        heatmap_maker(InN1_PV_data, "InN1_PV", marker = FALSE) %>% draw
dev.off()

InN9_PV_data <- filtered_results %>% filter(cluster_id == "InN9_PV") %>% arrange(logFC) %>%
                                select(any_of(c(metadata$Sample, "gene")))
write.csv(InN9_PV_data, file = paste0(base_path, "Finalized_outputs/14_InN9_PV_heatmap_source_data.csv"))
svg(file = paste0(base_path, "Finalized_outputs/14_Female_InN9_PV_heatmap.svg"), height = 6, width = 6)
        heatmap_maker(InN9_PV_data, "InN9_PV", marker = FALSE) %>% draw
dev.off()

InN8_ADARB2_data <- filtered_results %>% filter(cluster_id == "InN8_ADARB2") %>% arrange(logFC) %>%
                                select(any_of(c(metadata$Sample, "gene")))
write.csv(InN8_ADARB2_data, file = paste0(base_path, "Finalized_outputs/14_InN8_ADARB2_heatmap_source_data.csv"))
svg(file = paste0(base_path, "Finalized_outputs/14_Female_InN8_ADARB2_heatmap.svg"), height = 6, width = 6)
        heatmap_maker(InN8_ADARB2_data, "InN8_ADARB2", marker = FALSE) %>% draw
dev.off()

InN2_SST_data <- filtered_results %>% filter(cluster_id == "InN2_SST") %>% arrange(logFC) %>%
                                select(any_of(c(metadata$Sample, "gene")))
write.csv(InN2_SST_data, file = paste0(base_path, "Finalized_outputs/14_InN2_SST_heatmap_source_data.csv"))
svg(file = paste0(base_path, "Finalized_outputs/14_Female_InN2_SST_heatmap.svg"), height = 6, width = 8)
        heatmap_maker(InN2_SST_data, "InN2_SST", marker = FALSE) %>% draw
dev.off()

directory <- directories[["male_subtype"]]

filtered_results <- read_csv(file = paste0(base_path, directory, "/01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct.csv"))

#try this later: https://jokergoo.github.io/ComplexHeatmap-reference/book/genome-level-heatmap.html

ExN10_L46_data <- filtered_results %>% filter(cluster_id == "ExN10_L46") %>% arrange(logFC) %>%
                                select(any_of(c(metadata$Sample, "gene")))
write.csv(ExN10_L46_data, file = paste0(base_path, "Finalized_outputs/14_ExN10_heatmap_source_data.csv"))
svg(file = paste0(base_path, "Finalized_outputs/14_Male_ExN10_L46_heatmap.svg"), height = 16, width = 8)
        heatmap_maker(ExN10_L46_data, "ExN10_L46", marker = FALSE, h_dim = 0.1, f_size = 6, selected_labels = 20) %>% draw
dev.off()

Ast1_data <- filtered_results %>% filter(cluster_id == "Ast1") %>% arrange(logFC) %>%
                                select(any_of(c(metadata$Sample, "gene")))
write.csv(Ast1_data, file = paste0(base_path, "Finalized_outputs/14_Ast1_heatmap_source_data.csv"))
svg(file = paste0(base_path, "Finalized_outputs/14_Male_Ast1_heatmap.svg"), height = 12, width = 8)
        heatmap_maker(Ast1_data, "Ast1", marker = FALSE, h_dim = 0.2, f_size = 7) %>% draw
dev.off()

print(summary(warnings()))
sessionInfo()
