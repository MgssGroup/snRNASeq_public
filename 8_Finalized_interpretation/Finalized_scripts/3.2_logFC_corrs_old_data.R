library(tidyverse)
library(ggpmisc)
library(corrplot)
library(psych)

print("Set the base path and list of directories")
base_path <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/"

directories <- c("male_broad"=  "04_Male_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Male_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
                "male_subtype" = "Mar9_2022_updated_res/")

print("Read in old results files for previous male dataset analysis")

old_data_files <- list.files(paste0(base_path, "PreviousDiffExp_ForCorr"))

names(old_data_files) <- paste0("Prev_", str_sub(old_data_files, 1, -25))

print(old_data_files)

old_res_for_corr <- list()

for(file_name in names(old_data_files)) {
	print(file_name)
	results <- read_csv(file = paste0(base_path, "PreviousDiffExp_ForCorr/",  old_data_files[[file_name]]))
	colnames(results)[1] <- "gene"
	results %>% filter(p.adjust < 0.1) %>% select(gene) %>% print()
	results %>% mutate(score = Estimate * (-log10(`Pr(>|t|)`))) %>% select(gene, Estimate, score) -> res_for_corr
	colnames(res_for_corr)[2:3] <- paste(file_name, colnames(res_for_corr)[2:3], sep = ".")
	old_res_for_corr[[file_name]] <- res_for_corr

}

print("Read in new data files for the male dataset in the new analysis")
new_res_for_corr <- list()

for(directory_name in names(directories)) {
	print(directory_name)
        directory <- directories[[directory_name]]
	
	results <- read_csv(file = paste0(base_path, directory, "01_Table_raw_masterfile_cpm_Freq_DS_stats_greater3.csv")) %>% select(-c(1))

        print("Get the logFC and score for each gene per cluster")
        res_for_corr <- results %>% mutate(score = logFC * (-log10(p_val))) %>% group_by(cluster_id) %>% 
		#filter(grepl("Ast|OPC|ExN10|ExN3", cluster_id)) %>%
		select(gene, logFC, score, p_val) %>% group_split(.keep = FALSE)
        names(res_for_corr) <- results %>% group_by(cluster_id) %>% 
		#filter(grepl("Ast|OPC|ExN10|ExN3", cluster_id)) %>%
		group_keys() %>% select(cluster_id) %>% unlist()
	#print(res_for_corr)
        new_res_for_corr[[directory_name]] <- lapply(names(res_for_corr), function(x) {
		this_res <- res_for_corr[[x]]  
		setNames(this_res, nm = c(names(this_res)[1], paste(x, names(this_res)[2:4], sep = ".")))
        })
        names(new_res_for_corr[[directory_name]]) <- results %>% group_by(cluster_id) %>% 
		#filter(grepl("Ast|OPC|ExN10|ExN3", cluster_id)) %>%
                group_keys() %>% select(cluster_id) %>% unlist()

}

print("Function to plot fold change relationships for clusters of special interest")
scatter_plots_FCs <- function(df, combo, top = FALSE) {
	print(combo)
	print(top)
	p_val_col <- str_replace(combo[1], pattern = "logFC", replacement = "p_val") %>%
		str_replace(pattern = "score", replacement = "p_val")
	print(p_val_col)
	if(top == TRUE) {
		this_df <- df %>% select(any_of(c(combo, p_val_col))) %>% drop_na() %>% arrange(get(p_val_col))
		print(head(this_df))
                print(tail(this_df)) 
		this_df %>% slice_head(n = 1000) %>% select(any_of(combo)) -> this_df
                print(head(this_df))
                print(tail(this_df))
        } else {	
		this_df <- df %>% select(any_of(combo)) %>% drop_na
		print(head(this_df))
		print(tail(this_df))
	}
	print(summary(this_df[,combo[1]]))
	print(summary(this_df[,combo[2]]))	
        this_x = combo[1]
        this_y = combo[2]
        this_form = as.formula(paste0(this_y, "~", this_x))
        this_res <- lm(this_form, data = drop_na(this_df)) %>% summary()
        this_annot <- list(R_squared = signif(this_res$r.squared,3), p_val = signif(this_res$coefficients[2,4],3), df = this_res$df[2]) %>% as_tibble()
        if(grepl("Estimate", combo[2])) { label_x <- 1; label_y <- 1 }
	else { label_x <- 1.5; label_y <- 1.5} 
	p1 <- ((ggplot(data = drop_na(this_df), aes_string(x = this_x, y = this_y)) +
               geom_smooth(method = "lm", se=FALSE, color="black") +
               geom_point() +
               annotate("table", label = this_annot, x= label_x, y = label_y, size = 5) + #
               theme_classic(base_size = 22) + theme(axis.text = element_text(size = 22))) %>% ggrastr::rasterize())
	print(p1)
	write.csv(p1$data, file = paste0(base_path, "Finalized_outputs/3.2_source_data_", this_x, "_", this_y,".csv"))
}

print("Output correlations for broad clusters")
pdf(file = paste0(base_path, "Finalized_outputs/3.2_broad_new_old_corrs.pdf"), onefile = TRUE, height = 10, width = 10)
	corr_df <- reduce(c(new_res_for_corr$male_broad, old_res_for_corr), full_join, by = "gene") %>% column_to_rownames("gene")
	x_y_combos <- list(c("Ast.logFC", "Prev_Astros_3_AB.Estimate"), c("Ast.score", "Prev_Astros_3_AB.score"),
			c("Ast.logFC", "Prev_Astros_2.Estimate"), c("Ast.score", "Prev_Astros_2.score"),
			c("OPC.logFC", "Prev_OPCs_2.Estimate"), c("OPC.score", "Prev_OPCs_2.score"),
			c("OPC.logFC", "Prev_OPCs_1.Estimate"),  c("OPC.score", "Prev_OPCs_1.score"))
	for(combo in x_y_combos) {
		scatter_plots_FCs(corr_df, combo)
		scatter_plots_FCs(corr_df, combo, top = TRUE)
	}
	df_1 <- corr_df %>% select(contains("logFC"))
	df_2 <- corr_df %>% select(contains("Estimate"))
	corr_res <- corr.test(x = df_1, y = df_2, use = "pairwise", adjust = "BH")
	print(corr_res$n)
        corrplot(corr_res$r, p.mat = corr_res$p, tl.col = "black", insig = "label_sig", title = "Broad",
                mar = c(0,0,1,0), tl.cex = 0.7, pch.cex = 0.5) %>% print()
	pheatmap::pheatmap(corr_res$r)	
dev.off()

print("Output correlations for fine clusters")
pdf(file = paste0(base_path, "Finalized_outputs/3.2_cell_subtype_new_old_corrs.pdf"), onefile = TRUE)
	corr_df <- reduce(c(new_res_for_corr$male_subtype, old_res_for_corr), full_join, by = "gene") %>% column_to_rownames("gene")
	x_y_combos <- list(c("Ast1.logFC", "Prev_Astros_3_AB.Estimate"), c("Ast1.score", "Prev_Astros_3_AB.score"),
                        c("Ast1.logFC", "Prev_Astros_2.Estimate"), c("Ast1.score", "Prev_Astros_2.score"),
			c("Ast2.logFC", "Prev_Astros_3_AB.Estimate"), c("Ast2.score", "Prev_Astros_3_AB.score"),
                        c("Ast2.logFC", "Prev_Astros_2.Estimate"), c("Ast2.score", "Prev_Astros_2.score"),
                        c("OPC1.logFC", "Prev_OPCs_2.Estimate"), c("OPC1.score", "Prev_OPCs_2.score"),
                        c("OPC1.logFC", "Prev_OPCs_1.Estimate"),  c("OPC1.score", "Prev_OPCs_1.score"),
			c("OPC2.logFC", "Prev_OPCs_2.Estimate"), c("OPC2.score", "Prev_OPCs_2.score"),
                        c("OPC2.logFC", "Prev_OPCs_1.Estimate"),  c("OPC2.score", "Prev_OPCs_1.score"),
			c("ExN10_L46.logFC", "Prev_Ex_7_L4_6.Estimate"), c("ExN10_L46.score", "Prev_Ex_7_L4_6.score"),
                        c("ExN10_L46.logFC", "Prev_Ex_6_L4_6.Estimate"), c("ExN10_L46.score", "Prev_Ex_6_L4_6.score"),
                        c("ExN3_L46.logFC", "Prev_Ex_7_L4_6.Estimate"), c("ExN3_L46.score", "Prev_Ex_7_L4_6.score"),
                        c("ExN3_L46.logFC", "Prev_Ex_6_L4_6.Estimate"), c("ExN3_L46.score", "Prev_Ex_6_L4_6.score"))
        for(combo in x_y_combos) {
		scatter_plots_FCs(corr_df, combo)
		scatter_plots_FCs(corr_df, combo, top = TRUE)
	}
	df_1 <- corr_df %>% select(contains("logFC"))
        df_2 <- corr_df %>% select(contains("Estimate"))
        print(corr_res$n)
	corr_res <- corr.test(x = df_1, y = df_2, use = "pairwise", adjust = "BH")
        corrplot(corr_res$r, p.mat = corr_res$p, tl.col = "black", insig = "label_sig", title = "Broad",
                mar = c(0,0,1,0), tl.cex = 0.7, pch.cex = 0.5) %>% print()
	pheatmap::pheatmap(corr_res$r)
dev.off()

print(summary(warnings()))
sessionInfo()
