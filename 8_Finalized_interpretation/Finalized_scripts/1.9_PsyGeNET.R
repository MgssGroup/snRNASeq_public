library(tidyverse)
library(psygenet2r)

base_path <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/"

directories <- c("male_broad" = "04_Male_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Male_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
                "male_subtype" = "Mar9_2022_updated_res/",
                "female_broad" = "03_Female_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
                "female_subtype" = "03_Female_dataset_analysis/01_pseudobulk/02_cell_subtype/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_subtype/")

print("Make pie plots of numbers of DEGs per cluster after filtering")
dirs <- c("Down", "Up")
#pie_plots <- list()
summary_plot <- tibble()

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
                        pos = if_else(is.na(pos), Freq/2, pos)) %>%
                        mutate(cluster_id = paste(cluster_id, Freq, sep = " - "))
                #pie_plots[[paste0(direction, "|" , directory_name)]] <- ggplot(data, aes(x="", y=Freq, fill=cluster_id, label = )) +
                #        geom_bar(stat="identity", width=1, label = "") +
                #        coord_polar("y", start=0) +
                #        theme_void(base_size = 22) +
                #        ggtitle(paste0(direction, "|" , directory_name))
        }
	filtered_results %>% summarise(dataset = directory_name, 
				common_genes = n()-n_distinct(gene), 
				unique_genes = n_distinct(gene), 
				upregulated = n()-sum(logFC < 0), 
				downregulated = sum(logFC < 0)) %>% rbind(summary_plot) -> summary_plot
}

base_path_out <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/PsyGeNET/"

print(summary_plot)

summary_plot %>% select(1, 4:5) %>% pivot_longer(cols = where(is.numeric), names_to = "category", values_to = "percentage") -> unique_genes_plot
summary_plot %>% select(1:3) %>% pivot_longer(cols = where(is.numeric), names_to = "category", values_to = "percentage") -> up_down_genes_plot

#pdf(file = paste0(base_path_out, "1.9_pie_charts.pdf"), onefile = TRUE, width = 10, height = 5)
#	par(mar = c(2,2,2,2))
#        pie_plots
#dev.off()

pdf(file = paste0(base_path_out, "1.9_barplot_summaries.pdf"), onefile = TRUE, width = 10, height = 4)
	p1 <- ggplot(unique_genes_plot, aes(fill=category, x=percentage, y=dataset)) + 
	    geom_col(position="fill", stat="identity") + theme_classic(base_size = 22) + labs(x = NULL, y = NULL) +
		 theme(axis.text.x = element_text(size = 22, angle = 45, hjust = 1), legend.position="bottom")
	p1
	write.csv(p1$data, file = paste0(base_path_out, "1.9_barplot_summaries_1.csv"))
	p2 <- ggplot(up_down_genes_plot, aes(fill=category, x=percentage, y=dataset)) +
            geom_col(position="fill", stat="identity") + theme_classic(base_size = 22) + labs(x = NULL, y = NULL) +
		 theme(axis.text.x = element_text(size = 22, angle = 45, hjust = 1), legend.position="bottom")
	p2
	write.csv(p2$data, file = paste0(base_path_out, "1.9_barplot_summaries_2.csv"))
dev.off()

cluster_lists <- list("male_broad" = c("Ast", "OPC", "Oli"), "male_subtype" = c("Ast1", "ExN10_L46", "InN3_VIP", "ExN9_L23", "Ast2"), 
		 "female_broad" = c("Mic"), "female_subtype" = c("Mic1", "InN1_PV", "InN8_ADARB2", "InN9_PV", "InN2_SST"))

for(directory_name in names(directories)) {
	pdf(file = paste0(base_path_out, "1.9_", directory_name,"_psygenet_GDA.pdf"), onefile = TRUE, width = 16, height = 10)
	results_list <- list()
	tables_list <- list()
	directory <- directories[[directory_name]]
	results <- read_csv(file = paste0(base_path, directory, "01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct.csv")) %>% select(-c(1))
	genesOfInterest <- results$gene %>% unique()
	print(paste("Running PsyGeNET overall for", directory_name))
	results_list[["all"]] <- psygenetGene(
		gene     = genesOfInterest, 
		database = "ALL",
		verbose  = TRUE)
	tables_list[["all"]] <- enrichedPD( genesOfInterest, database = "ALL")
	print(paste("Running PsyGeNET by cluster for", directory_name))
	for(cluster in cluster_lists[[directory_name]]) {
		genesOfInterest <- results %>% filter(cluster_id == cluster) %>% select(gene) %>% unlist()
		results_list[[cluster]] <- psygenetGene(
	                gene     = genesOfInterest,
	                database = "ALL",
	                verbose  = TRUE) 
		tables_list[[cluster]] <- enrichedPD(genesOfInterest, database = "ALL")
	}
	print("Plotting PsyGeNET results")
	lapply(names(results_list), function(x) {
		res_p <- plot(results_list[[x]], type = "GDA heatmap") + theme_classic(base_size = 22) + ggtitle(paste(directory_name, x)) +
                        theme(axis.text.y = element_text(size = 22), axis.text.x = element_text(size = 22, angle = 45, hjust = 1))
		print(res_p)
		write.csv(res_p$data, file = paste0(base_path_out, "1.9_source_data", directory_name, "_", x, "_psygenet_GDA_heatmap.csv"))
		res_p <- plot(results_list[[x]], type = "GDCA heatmap") + ggtitle(paste(directory_name, x)) + theme_classic(base_size = 22) +
                        theme(axis.text.y = element_text(size = 22), axis.text.x = element_text(size = 22, angle = 45, hjust = 1))
		print(res_p)
		res_p <- geneAttrPlot(results_list[[x]], type = "evidence index") + ggtitle(paste(directory_name, x)) + theme_classic(base_size = 22) +
                        theme(axis.text.y = element_text(size = 22), axis.text.x = element_text(size = 22, angle = 45, hjust = 1))
		print(res_p)
		write.csv(res_p$data, file = paste0(base_path_out, "1.9_source_data", directory_name, "_", x, "_psygenet_geenAttr.csv"))
	})
	dev.off()
	write.csv(bind_rows(tables_list, .id = "cluster_id"), file = paste0(base_path_out, "1.9_", directory_name, "_psygenet_table.csv"))
}

print(summary(warnings()))
sessionInfo()
