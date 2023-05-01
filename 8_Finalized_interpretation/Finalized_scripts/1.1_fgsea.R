library(data.table)
library(fgsea)
library(msigdbr)
library(tidyverse)
library(ggrepel)

base_path <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/"

directories <- c(male_broad = "04_Male_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Male_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
                male_subtype = "Mar9_2022_updated_res/",
                female_broad =  "03_Female_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
                female_subtype = "03_Female_dataset_analysis/01_pseudobulk/02_cell_subtype/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_subtype/")

print("Make pie plots of numbers of DEGs per cluster after filtering")
dirs <- c("Down", "Up")
pie_plots <- list()

for(directory_name in names(directories)) {
	directory <- directories[[directory_name]]
	filtered_results <- read_csv(file = paste0(base_path, directory, "01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct.csv")) %>% select(-c(1)) 
	for(direction in c("up", "down", "all", "combined")) {
                print(direction)
                data <- switch(direction,
                        "up" = filtered_results %>% filter(logFC > 0) %>% select(cluster_id) %>% unlist() %>% table() %>% as.data.frame(),
                        "down" = , filtered_results %>% filter(logFC < 0) %>% select(cluster_id) %>% unlist() %>% table() %>% as.data.frame(),
                        "all" = filtered_results %>% select(cluster_id) %>% unlist() %>% table() %>% as.data.frame(),
                        "combined" = filtered_results %>% mutate(dir = map_chr(logFC, ~dirs[sign(.x)/2+1.5])) %>%
                                mutate(cluster_id = paste0(dir, "_", cluster_id)) %>%
                                select(cluster_id) %>% unlist() %>% table() %>% as.data.frame())
                print(dim(data))
                colnames(data) <- c("cluster_id", "Freq")
		data <- data %>%
			mutate(csum = rev(cumsum(rev(Freq))), 
			pos = Freq/2 + lead(csum, 1),
			pos = if_else(is.na(pos), Freq/2, pos))
                pie_plots[[paste0(direction, "|" , directory_name)]] <- ggplot(data, aes(x="", y=Freq, fill=cluster_id)) +
                        geom_bar(stat="identity", width=1, color = "white") +
                        coord_polar("y", start=0) +
			theme_void() + 
			geom_label_repel(aes(y = pos, label = Freq), color = "black", size=6, show.legend = FALSE, nudge_x = 1) +
                        ggtitle(paste0(direction, "|" , directory_name))
        }
	
}

pdf(file = paste0(base_path, "Finalized_outputs/1_pie_charts.pdf"), onefile = TRUE)
        pie_plots
dev.off()

rm(pie_plots, filtered_results)

print("Run FGSEA analysis")

print("Define function for formatting misgdbr gene sets")
format_gene_set <- function(x, y) {
	#print(y$gs_name)
	formatted_gene_set_list <- list()
	formatted_gene_set_list[[y$gs_name]] <- x$gene_symbol
	formatted_gene_set_list
}

print("Get and format msigdbr gene sets")
subcategory_category <- list("MIR:MIRDB" = "C3", "TFT:GTRD" = "C3", "GO:BP" = "C5", "GO:MF" = "C5", "GO:CC" = "C5",  "HPO" = "C5", "CP:REACTOME" = "C2")

gene_sets <- lapply(names(subcategory_category), function(x) {
		 print(paste(x, subcategory_category[[x]]))	
		 msigdbr(species = "Homo sapiens", category = subcategory_category[[x]], subcategory = x) %>% 
			group_by(gs_name) %>% 
			group_map(~format_gene_set(.x, .y)) %>%
			unlist(recursive = FALSE)
})

names(gene_sets) <- names(subcategory_category)

print("Function for prepping genes stats lists per cluster")
prep_gene_stats <- function(x, y){
	#print(cluster_id)
	prepped_gene_stats <- x$score
	names(prepped_gene_stats) <- x$gene
	prepped_gene_stats_list <- list()
	prepped_gene_stats_list[[y$cluster_id]] <- prepped_gene_stats
	prepped_gene_stats_list	
}

print("Function to plot top enriched pathways for each cluster in each analysis.")
fgsea_results_plotting <- function(fgseaRes, pathways, ranks, file_name) {
	print("Write results to file")
        fwrite(fgseaRes, file=paste0(file_name, ".txt"), sep="\t", sep2=c("", " ", ""))

	print("Find main, collapsed pathways")
	collapsedPathways <- collapsePathways(fgseaRes[order(padj)][padj < 0.01], 
                                      pathways, ranks)
	mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
                         order(-NES), pathway]

	print("Find up and downregulated pathways")
	topPathwaysUp <- fgseaRes[ES > 0 & padj < 0.05][head(order(padj)), pathway]
        topPathwaysDown <- fgseaRes[ES < 0 & padj < 0.05][head(order(-padj)), pathway]
        topPathways <- c(topPathwaysUp, topPathwaysDown)
	
	print("Save plots to files")
	if(
		#length(mainPathways) > 0 | 
		length(topPathways) > 0) {
		pdf(file = paste0(file_name, ".pdf"), height = 10, width = 14, onefile = TRUE)
		#if(length(mainPathways) > 0) {
		#	plotGseaTable(pathways[mainPathways], ranks, fgseaRes, 
		#		gseaParam = 0.5) 
		#	while(!par('page')) plot.new()}	
		#if(length(topPathways) > 0) {
			plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
				gseaParam=0.5)
		#}
		dev.off()
	}
	print("Find up and downregulated pathways with more lenient threshold")
        topPathwaysUp <- fgseaRes[ES > 0 & padj < 0.1][order(padj), ]
        topPathwaysDown <- fgseaRes[ES < 0 & padj < 0.1][order(padj), ]
        topPathways <- rbind(topPathwaysUp, rev(topPathwaysDown))

        print("Save lenient top pathways  to file")
        if(length(topPathways[, pathway]) > 0) fwrite(topPathways, file=paste0(file_name, "_significant.txt"), sep="\t", sep2=c("", " ", ""))
	NULL
}

print("Run fgsea for each analysis, for each cluster, for each category of gene sets")
for(directory_name in names(directories)) {
	print(directory_name)
        directory <- directories[[directory_name]] 
	results <- read_csv(file = paste0(base_path, directory, "01_Table_raw_masterfile_cpm_Freq_DS_stats_greater3.csv")) %>% select(-c(1))
	print("Add ranking score")
	results %>% mutate(score = logFC * (-log10(p_val))) %>% 
		group_by(cluster_id) %>% select(gene, score) %>% 
                arrange(score, .by_group = TRUE) %>% group_map(~prep_gene_stats(.x, .y)) %>% unlist(recursive = FALSE) -> gene_stats
	for(gene_set in names(gene_sets)) {
		print(gene_set)
		pathways <- gene_sets[[gene_set]]
		res_fgsea  <- lapply(gene_stats, function(x) {
			set.seed(42)
			fgsea(pathways = pathways, 
				stats = x,
				eps      = 0.0,
				minSize  = 15,
				maxSize  = 1000, 
				nproc = 1)
		})
		dir_name <- paste0(directory_name, "_", gene_set)
		dir.create(paste0(base_path, "Finalized_outputs/FGSEA_results/",dir_name))
		temp <- lapply(names(res_fgsea), function(x) {
			ranks <- gene_stats[[x]]
			file_name <- paste0(base_path, "Finalized_outputs/FGSEA_results/", dir_name, "/", directory_name, "_", x, "_", gene_set)
			fgsea_results_plotting(res_fgsea[[x]], pathways, ranks, file_name)
		})
		rm(temp)	
	}
}

print(summary(warnings()))
sessionInfo()
