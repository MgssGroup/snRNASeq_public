library(tidyverse)
library(MASS)

print("Load original filtered files and count numbers of DEGs, clusters with DEGs, and overlaps")

original_directory <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/"

male_subtype <- read_csv(paste0(original_directory, 
	"Mar9_2022_updated_res/01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct.csv")) %>%
	 unite(ClusterGene, c("cluster_id", "gene"), remove = FALSE)
male_broad <- read_csv(paste0(original_directory, 
	"04_Male_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Male_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct.csv")) %>%
	 unite(ClusterGene, c("cluster_id", "gene"), remove = FALSE)
female_subtype <- read_csv(paste0(original_directory, 
	"03_Female_dataset_analysis/01_pseudobulk/02_cell_subtype/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_subtype/01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct.csv")) %>%
	 unite(ClusterGene, c("cluster_id", "gene"), remove = FALSE) 
female_broad <- read_csv(paste0(original_directory, 
	"03_Female_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct.csv")) %>%
	 unite(ClusterGene, c("cluster_id", "gene"), remove = FALSE) 

male_subtype %>% summarise(num_DEGs = n_distinct(gene), 
	num_cluster = n_distinct(cluster_id),
	num_cluster_gene = n_distinct(ClusterGene)) -> male_sub_nums
female_subtype %>% summarise(num_DEGs = n_distinct(gene), 
	num_cluster = n_distinct(cluster_id), 
	num_cluster_gene = n_distinct(ClusterGene)) -> female_sub_nums
female_broad %>% summarise(num_DEGs = n_distinct(gene), 
	num_cluster = n_distinct(cluster_id),
	num_cluster_gene = n_distinct(ClusterGene)) -> female_broad_nums
male_broad %>% summarise(num_DEGs = n_distinct(gene), 
	num_cluster = n_distinct(cluster_id),
	num_cluster_gene = n_distinct(ClusterGene)) -> male_broad_nums
overlap_broad_nums <- tibble(
	num_DEGs = length(intersect(male_broad$gene, female_broad$gene)), 
	num_cluster = length(intersect(male_broad$cluster_id, female_broad$cluster_id)),
	num_cluster_gene = length(intersect(male_broad$ClusterGene, female_broad$ClusterGene)))
overlap_subtypes_nums <- tibble(
	num_DEGs = length(intersect(male_subtype$gene, female_subtype$gene)), 
	num_cluster = length(intersect(male_subtype$cluster_id, female_subtype$cluster_id)),
	num_cluster_gene = length(intersect(male_subtype$ClusterGene, female_subtype$ClusterGene)))
original_numbers <- bind_rows(female_broad_nums, female_sub_nums, male_broad_nums, male_sub_nums, overlap_broad_nums, overlap_subtypes_nums)
original_numbers$type <- c("fb_nums", "fs_nums", "mb_nums", "ms_nums", "ob_nums", "os_nums")

rm(list = ls(pattern = "*nums"))
rm(list = ls(pattern = "male*"))
rm(list = ls(pattern = "female*"))

permutation_directory <- "/home/malosree/projects/def-cnagy/malosree/Permutation_outputs/"

subtype_nums <- list()
broad_nums <- list()

for(i in 1:100) {
	ms <- read_csv(paste0(permutation_directory, "1_male_subtype/1_male_subtype/", i, "_unfiltered_male_subtype.csv"))
	print(paste("Filter male subtype DEG table", i))
        ms %>% mutate(NumNonZero = rowSums(across(ends_with(".cpm"), ~(.x>0)), na.rm = TRUE)) %>%
                        mutate(Greater3 = NumNonZero > 2) -> ms_filtered
        ms_filtered %>% filter(p_adj.loc < 0.05, abs(logFC) > log2(1.1), Greater3 == TRUE) %>%
		unite(ClusterGene, c("cluster_id", "gene"), remove = FALSE)-> ms_filtered
	ms_filtered %>% summarise(num_DEGs = n_distinct(gene), 
				num_cluster = n_distinct(cluster_id), 
				num_cluster_gene = n_distinct(ClusterGene)) -> ms_nums
	rm(ms)	

	fs <- read_csv(paste0(permutation_directory, "3_female_subtype/",i, "_unfiltered_female_subtype.csv"))
	print(paste("Filter female subtype DEG table", i))
        fs %>% mutate(NumNonZero = rowSums(across(ends_with(".cpm"), ~(.x>0)), na.rm = TRUE)) %>%
                        mutate(Greater3 = NumNonZero > 2) -> fs_filtered
        fs_filtered %>% filter(p_adj.loc < 0.05, abs(logFC) > log2(1.1), Greater3 == TRUE) %>%
		unite(ClusterGene, c("cluster_id", "gene"), remove = FALSE)-> fs_filtered
        rm(fs)
	fs_filtered %>% summarise(num_DEGs = n_distinct(gene), 
				num_cluster = n_distinct(cluster_id),
				num_cluster_gene = n_distinct(ClusterGene)) -> fs_nums

	os_nums <- tibble(num_DEGs = length(intersect(ms_filtered$gene, fs_filtered$gene)), 
			num_cluster = length(intersect(ms_filtered$cluster_id, fs_filtered$cluster_id)),
			num_cluster_gene = length(intersect(ms_filtered$ClusterGene, fs_filtered$ClusterGene)))
	
	subtype_nums[[i]] <- bind_rows(fs_nums, ms_nums, os_nums)
	subtype_nums[[i]]$type <- c("fs_nums", "ms_nums", "os_nums")

	mb <- read_csv(paste0(permutation_directory, "2_male_broad/", i, "_unfiltered_male_broad.csv"))
        print(paste("Filter male broad DEG table", i))
        mb %>% mutate(NumNonZero = rowSums(across(ends_with(".cpm"), ~(.x>0)), na.rm = TRUE)) %>%
                        mutate(Greater3 = NumNonZero > 2) -> mb_filtered
        mb_filtered %>% filter(p_adj.loc < 0.05, abs(logFC) > log2(1.1), Greater3 == TRUE) %>%
			unite(ClusterGene, c("cluster_id", "gene"), remove = FALSE)-> mb_filtered
        mb_filtered %>% summarise(num_DEGs = n_distinct(gene), 
				num_cluster = n_distinct(cluster_id), 
				num_cluster_gene = n_distinct(ClusterGene)) -> mb_nums
        rm(mb)

        fb <- read_csv(paste0(permutation_directory, "4_female_broad/",i, "_unfiltered_female_broad.csv"))
        print(paste("Filter female broad DEG table", i))
        fb %>% mutate(NumNonZero = rowSums(across(ends_with(".cpm"), ~(.x>0)), na.rm = TRUE)) %>%
                        mutate(Greater3 = NumNonZero > 2) -> fb_filtered
        fb_filtered %>% filter(p_adj.loc < 0.05, abs(logFC) > log2(1.1), Greater3 == TRUE) %>% 
			unite(ClusterGene, c("cluster_id", "gene"), remove = FALSE)-> fb_filtered
        rm(fb)
        fb_filtered %>% summarise(num_DEGs = n_distinct(gene), 
				num_cluster = n_distinct(cluster_id), 
				num_cluster_gene = n_distinct(ClusterGene)) -> fb_nums

        ob_nums <- tibble(num_DEGs = length(intersect(mb_filtered$gene, fb_filtered$gene)),
                        num_cluster = length(intersect(mb_filtered$cluster_id, fb_filtered$cluster_id)), 
			num_cluster_gene = length(intersect(mb_filtered$ClusterGene,  fb_filtered$ClusterGene)))

        broad_nums[[i]] <- bind_rows(fb_nums, mb_nums, ob_nums)
        broad_nums[[i]]$type <- c("fb_nums", "mb_nums", "ob_nums")

}

broad_nums %>% bind_rows(.id = "iteration") -> broad_nums
subtype_nums %>% bind_rows(.id = "iteration") -> subtype_nums

print("Make distribution plots")


labels <- c(fb_nums = "female broad", 
		mb_nums = "male broad", 
		fs_nums = "female cluster",
		ms_nums = "male cluster", 
		ob_nums = "overlapping broad", 
		os_nums = "overlapping cluster")

fractions_below <- tibble()

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Finalized_outputs/5_distributions_individual_overlaps.pdf", 
	onefile = TRUE)
	for(type in unique(broad_nums$type)) {
		print(type)
		if(type != "ob_nums") {h_val <- "h = 50"} else {h_val <- ""}
		par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
		eval(parse(text = paste0('truehist(broad_nums$num_DEGs[broad_nums$type == type], xlab = paste("DEGs", labels[type]), prob = FALSE,', h_val,')'))) 
		abline(v = original_numbers$num_DEGs[original_numbers$type == type], col = "red")
		fractions_below <- bind_rows(fractions_below, tibble(this_type = type,
                        fraction = cume_dist(
                                c(original_numbers$num_DEGs[original_numbers$type == type],
                                        broad_nums$num_DEGs[broad_nums$type == type]))[1]))
		write.csv(broad_nums$num_DEGs[broad_nums$type == type], file = paste0(
			"/home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Finalized_outputs/5_source_data",
			type, "permuted_num_DEGs.csv"))		
		par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
		truehist(broad_nums$num_cluster[broad_nums$type == type], xlab = paste("Clusters", labels[type]), prob = FALSE, h = 1) 
		abline(v = original_numbers$num_cluster[original_numbers$type == type], col = "red")
		par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
                truehist(broad_nums$num_cluster_gene[broad_nums$type == type], 
			xlab = paste("Cluster gene combinations", labels[type]), prob = FALSE)
                abline(v = original_numbers$num_cluster_gene[original_numbers$type == type], col = "red")
	}
	for(type in unique(subtype_nums$type)) {
		print(type)
		if(type != "os_nums") {h_val <- "h = 50"} else { h_val <- ""}
		par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
		eval(parse(text = paste0('truehist(subtype_nums$num_DEGs[subtype_nums$type == type], xlab = paste("DEGs", labels[type]), prob = FALSE,', h_val, ')'))) 
		abline(v = original_numbers$num_DEGs[original_numbers$type == type], col = "red")
		par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
		truehist(subtype_nums$num_cluster[subtype_nums$type == type], xlab = paste("Clusters", labels[type]), prob = FALSE, h = 1) 
		abline(v = original_numbers$num_cluster[original_numbers$type == type], col = "red")
		par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
                truehist(subtype_nums$num_cluster_gene[subtype_nums$type == type], 
			xlab = paste("Cluster gene combinations", labels[type]), prob = FALSE)
                abline(v = original_numbers$num_cluster_gene[original_numbers$type == type], col = "red")
		fractions_below <- bind_rows(fractions_below, tibble(this_type = type,
                        fraction = cume_dist(
                                c(original_numbers$num_DEGs[original_numbers$type == type],
                                        subtype_nums$num_DEGs[subtype_nums$type == type]))[1]))
		write.csv(subtype_nums$num_DEGs[subtype_nums$type == type], file = paste0(
                        "/home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Finalized_outputs/5_source_data",
                        type, "permuted_num_DEGs.csv"))

        }
dev.off()
write.csv(fractions_below, file = 
	"/home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Finalized_outputs/5_DEG_fractions_permutations.csv")
print(summary(warnings()))
sessionInfo()
