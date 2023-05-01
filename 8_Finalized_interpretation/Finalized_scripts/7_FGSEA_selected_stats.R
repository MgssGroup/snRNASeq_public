library(Seurat)
library(msigdbr)
library(tidyverse)
library(AUCell)
library(rstatix)
library(SummarizedExperiment)
library(fgsea)
library(data.table)
library(edgeR)

base_path <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/"

directories <- c(
		#male_broad = "04_Male_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Male_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
                #male_subtype = "Mar9_2022_updated_res/",
                #female_broad =  "03_Female_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
                female_subtype = "03_Female_dataset_analysis/01_pseudobulk/02_cell_subtype/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_subtype/")

cells_to_use <- list(female_subtype = c(sex = "Female", class = "Cluster"))

clusters_of_interest <- list(female_subtype = c("Mic1", "InN1_PV", "InN9_PV"))

harmonized_object <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

print("Define function for formatting misgdbr gene sets")
format_gene_set <- function(x, y) {
        #print(y$gs_name)
        formatted_gene_set_list <- list()
        formatted_gene_set_list[[y$gs_name]] <- x$gene_symbol
        formatted_gene_set_list
}

print("Get and format msigdbr gene sets")
subcategory_category <- list(
				#"MIR:MIRDB" = "C3", "TFT:GTRD" = "C3", "GO:BP" = "C5", "GO:MF" = "C5", "GO:CC" = "C5",  "HPO" = "C5", 
				"CP:REACTOME" = "C2")

gene_sets <- lapply(names(subcategory_category), function(x) {
                 print(paste(x, subcategory_category[[x]]))
                 msigdbr(species = "Homo sapiens", category = subcategory_category[[x]], subcategory = x) %>%
                        group_by(gs_name) %>%
                        group_map(~format_gene_set(.x, .y)) %>%
                        unlist(recursive = FALSE)
})

names(gene_sets) <- names(subcategory_category)

print("Functon for wilcoxon tests for subject geneset scores between cases and controls")
wilcoxon_genesets <- function(geneset_values, conditions){
        myData <- data.frame(GeneSetScore = geneset_values, Condition = as.factor(conditions))
        res_test <- wilcox_test(formula = GeneSetScore ~ Condition, data = myData, alternative = "two.sided", paired = FALSE, detailed = TRUE)
}

#Add module scores and do stats 
for(dir_name in names(directories)) {
	for(set_name  in names(gene_sets)) {
		for(cell_type in clusters_of_interest[[dir_name]]) {
			combination_name <- paste0(dir_name, "_", cell_type, "_", set_name)
                        print(combination_name)
			pathways <- read_tsv(file = paste0(base_path, "Finalized_outputs/FGSEA_results/", dir_name, "_",
                                        set_name, "/", combination_name, "_significant.txt"))				
			#Add, plot, stats for Seurat module scores
			obj <- subset(harmonized_object, Sex == cells_to_use[[dir_name]]["sex"] &
                                                harmonized_object@meta.data[cells_to_use[[dir_name]]["class"]] == cell_type)
			obj <- AddModuleScore(obj, features = gene_sets[[set_name]][unlist(pathways$pathway)])
			for(i in 1:length(pathways$pathway)) {
				obj@meta.data[, pathways$pathway[i]] <- obj@meta.data[,paste0("Cluster", i)]
				obj@meta.data[,paste0("Cluster", i)] <- NULL
			}
			seurat_to_plot <- obj@meta.data %>% group_by(Sample, Condition) %>% select(pathways$pathway) %>%
				summarise_all(mean) %>% ungroup()
			lapply(select(seurat_to_plot, where(is.numeric)), function(x) {
                                        wilcoxon_genesets(x, seurat_to_plot$Condition)
                                }) %>% bind_rows(.id = "Cluster")-> results_wilcoxon_seurat
			results_wilcoxon_seurat %>% mutate(p_adjust = p.adjust(p, method = "BH")) -> results_wilcoxon_seurat
			pdf(file = paste0(base_path, "Finalized_outputs/Geneset_stats/7_Seurat_scores_", combination_name, ".pdf"))
			for(column_name in names(select_if(seurat_to_plot, is.numeric))) {
                                (ggplot(seurat_to_plot, aes(y = seurat_to_plot[[column_name]], x = Condition, fill = Condition, label = Sample)) + 
					geom_boxplot() + theme_classic() +
                                        geom_label() + ylab(column_name)) %>%
                                        ggrastr::rasterize() %>% print
                        }
			dev.off()	 
			write.csv(results_wilcoxon_seurat, file = paste0(base_path, 
					"Finalized_outputs/Geneset_stats/7_AUCell_wilcoxon_seurat", combination_name, ".csv"))
		}
	}
}

#Collapse pathways and run camera

for(dir_name in names(directories)) {
        for(set_name  in names(gene_sets)) {
                for(cell_type in clusters_of_interest[[dir_name]]) {
                        combination_name <- paste0(dir_name, "_", cell_type, "_", set_name)
                        print(combination_name)
                        fGseaRes <- read_tsv(file = paste0(base_path, "Finalized_outputs/FGSEA_results/", dir_name, "_",
                                        set_name, "/", combination_name, "_significant.txt")) %>% as.data.frame() %>% setDT()			
			pathways <- gene_sets[[set_name]]
			results <- read_csv(file = paste0(base_path, directories[[dir_name]], "01_Table_raw_masterfile_cpm_Freq_DS_stats_greater3.csv")) %>% select(-c(1))
			print("Add ranking score")
			results %>% filter(cluster_id == cell_type) %>% mutate(score = logFC * (-log10(p_val))) %>%
				select(gene, score) %>% arrange(score) %>% deframe()-> gene_stats
			collapsed_pathways <- collapsePathways(fGseaRes, pathways, gene_stats, pval.threshold= 0.01) 
			dge_res <-  readRDS(paste0(base_path, directories[[dir_name]], "res.rds"))
			dge_res <- dge_res$data[[cell_type]]
			ids2indices(pathways[collapsed_pathways$mainPathways], rownames(dge_res)) -> pathway_indices
			camera_res_collapsed <- camera.DGEList(dge_res, pathway_indices)
			write.csv(collapsed_pathways$parentPathways, file = paste0(base_path,
                                        "Finalized_outputs/Geneset_stats/7_collapsed_pathways_", combination_name, ".csv"))
			write.csv(camera_res_collapsed, file = paste0(base_path,
                                        "Finalized_outputs/Geneset_stats/7_camera_test_collapsed_pathways_", combination_name, ".csv"))
			if(set_name == "CP:REACTOME") {
				reactome_pathways <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
				reactome_pathways %>% filter(gs_name %in% collapsed_pathways$mainPathways) %>% 
					select(gs_exact_source,gs_name) %>% distinct -> reactome_ids
				write.csv(reactome_ids, file = paste0(base_path,
                                        "Finalized_outputs/Geneset_stats/7_collapsed_pathways_", combination_name, "_reactome_ids.csv"))
			}
                }
        }
}

print(summary(warnings()))
sessionInfo()
