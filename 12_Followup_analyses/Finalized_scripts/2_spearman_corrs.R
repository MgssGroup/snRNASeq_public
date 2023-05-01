library(tidyverse)

#Broad#####

setwd("/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/")

res_male_raw <-
  read_csv(file = 
	"04_Male_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Male_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/01_Table_raw_masterfile_cpm_Freq_DS_stats_greater3.csv") %>% 
	as.data.frame()
res_male_raw$ID <-
  paste0(res_male_raw$gene, "_", res_male_raw$cluster_id)
res_male <- select(res_male_raw, ID, p_val, logFC)	
res_male$Score  <-  -log10(res_male$p_val)
res_male$Score <- res_male$Score * (res_male$logFC)

gene_list1 <-
  data.frame(Genes = res_male$ID, DDE = res_male$Score)

res_female_raw <-
  read_csv(
    "03_Female_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/01_Table_raw_masterfile_cpm_Freq_DS_stats_greater3.csv") %>% 
	as.data.frame()
res_female_raw$ID <-
  paste0(res_female_raw$gene, "_", res_female_raw$cluster_id)
res_female <- select(res_female_raw, ID, p_val, logFC)
res_female$Score  <-  -log10(res_female$p_val)
res_female$Score <- res_female$Score * (res_female$logFC)

gene_list2 <-
  data.frame(Genes = res_female$ID, DDE = res_female$Score)

df <- inner_join(gene_list1, gene_list2,  by = c("Genes"), suffix = c(".male", ".female"))
gene_list1 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.male)
gene_list2 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.female)


#gene_list1 <- sample_n(gene_list1, 20000)

corr_results <- data.frame()

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_outputs/2_rankplots_broad.pdf", onefile = TRUE)

clusterID <- intersect(unique(res_female_raw$cluster_id), unique(res_male_raw$cluster_id)) %>% sort()
for (i in 1:length(clusterID)) {
  #i = 2
  print(paste0("Running clusterID = ", clusterID[i]))
  #i =1
  gene_list1_cluster <- gene_list1[grep(paste0("_", clusterID[i]),gene_list1$Genes),]
  gene_list2_cluster <- gene_list2[grep(paste0("_", clusterID[i]),gene_list2$Genes),]
  gene_list2_cluster <- gene_list2_cluster[match(gene_list1_cluster$Genes, gene_list2_cluster$Genes), ]
  
  gene_list_1_2_cluster <- inner_join(gene_list1_cluster, gene_list2_cluster, by = c("Genes"), suffix = c(".male", ".female"))

  cor_res <- cor.test(gene_list_1_2_cluster$DDE.female, gene_list_1_2_cluster$DDE.male, method = "spearman")	
  cor_res <- c(clusterID[i], cor_res[c("estimate", "statistic", "p.value", "alternative", "method")])	
  names(cor_res) <- c("cluster", "estimate", "statistic", "p.value", "alternative", "method")

  corr_results <- rbind(corr_results, cor_res)	

  (ggplot(gene_list_1_2_cluster, aes(x = DDE.female, y = DDE.male)) + geom_point() + theme_classic(base_size=20) + ggtitle(clusterID[i])) %>% 
	ggrastr::rasterize() %>% print

  gene_list_1_2_cluster$rank_male <- rank(-gene_list_1_2_cluster$DDE.male)
  gene_list_1_2_cluster$rank_female <- rank(-gene_list_1_2_cluster$DDE.female)
	
  (gene_list_1_2_cluster$DDE.male > 0)	%>% sum -> x_line
  (gene_list_1_2_cluster$DDE.female > 0)  %>% sum -> y_line

  #(ggplot(gene_list_1_2_cluster, aes(x = rank_female, y = rank_male)) + geom_point() + theme_classic(base_size=20) + ggtitle(clusterID[i])) %>%
  #		ggrastr::rasterize() %>% print
  #(ggplot(gene_list_1_2_cluster, aes(x = rank_female, y = rank_male)) + geom_hex(binwidth = RRHO2:::defaultStepSize(gene_list1_cluster, gene_list2_cluster)) +
  #              scale_fill_continuous(type = "viridis") + theme_classic(base_size=20) +
  #              geom_vline(xintercept = y_line) + geom_hline(yintercept = x_line) + ggtitle(clusterID[i])) %>%  ggrastr::rasterize() %>% print
  (ggplot(gene_list_1_2_cluster, aes(x = rank_female, y = rank_male)) + geom_density_2d_filled(show.legend = FALSE) + theme_classic(base_size=20) + 
		geom_vline(xintercept = y_line) + geom_hline(yintercept = x_line) + ggtitle(clusterID[i])) %>%  ggrastr::rasterize() %>% print
	
}

dev.off()

names(corr_results) <- c("cluster", "estimate", "statistic", "p.value", "alternative", "method")

write.csv(corr_results, file = "/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_outputs/2_spearman_broad.csv")

#cell_subtype####

res_male_raw <-
  read_csv(file = "Mar9_2022_updated_res/01_Table_raw_masterfile_cpm_Freq_DS_stats_greater3.csv") %>% as.data.frame()
res_male_raw$ID <-
  paste0(res_male_raw$gene, "_", res_male_raw$cluster_id)
res_male <- select(res_male_raw, ID, p_val, logFC)
res_male$Score  <- -log10(res_male$p_val)
res_male$Score <- res_male$Score * (res_male$logFC)

gene_list1 <-
  data.frame(Genes = res_male$ID, DDE = res_male$Score)

res_female_raw <-
  read_csv(
    "03_Female_dataset_analysis/01_pseudobulk/02_cell_subtype/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_subtype/01_Table_raw_masterfile_cpm_Freq_DS_stats_greater3.csv") %>% 
	as.data.frame()
res_female_raw$ID <-
  paste0(res_female_raw$gene, "_", res_female_raw$cluster_id)
res_female <- select(res_female_raw, ID, p_val, logFC)
res_female$Score  <-  -log10(res_female$p_val)
res_female$Score <- res_female$Score * (res_female$logFC)

gene_list2 <-
  data.frame(Genes = res_female$ID, DDE = res_female$Score)

df <- inner_join(gene_list1, gene_list2,  by = c("Genes"), suffix = c(".male", ".female"))
gene_list1 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.male)
gene_list2 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.female)

corr_results <- data.frame()

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_outputs/2_rankplots_cluster.pdf", onefile = TRUE)

clusterID <- intersect(unique(res_female_raw$cluster_id), unique(res_male_raw$cluster_id)) %>% sort()
for (i in 1:length(clusterID)) {
  #i = 2
  print(paste0("Running clusterID = ", clusterID[i]))
  #i =1
  gene_list1_cluster <- gene_list1[grep(paste0("_", clusterID[i]),gene_list1$Genes),]
  gene_list2_cluster <- gene_list2[grep(paste0("_", clusterID[i]),gene_list2$Genes),]
  gene_list2_cluster <- gene_list2_cluster[match(gene_list1_cluster$Genes, gene_list2_cluster$Genes), ]

  gene_list_1_2_cluster <- inner_join(gene_list1_cluster, gene_list2_cluster, by = c("Genes"), suffix = c(".male", ".female"))

  cor_res <- cor.test(gene_list_1_2_cluster$DDE.female, gene_list_1_2_cluster$DDE.male, method = "spearman")
  cor_res <- c(clusterID[i], cor_res[c("estimate", "statistic", "p.value", "alternative", "method")])
  names(cor_res) <- c("cluster", "estimate", "statistic", "p.value", "alternative", "method")

  corr_results <- rbind(corr_results, cor_res)

  (ggplot(gene_list_1_2_cluster, aes(x = DDE.female, y = DDE.male)) + geom_point() + theme_classic(base_size=20) + ggtitle(clusterID[i])) %>% 
	ggrastr::rasterize() %>% print

  gene_list_1_2_cluster$rank_male <- rank(-gene_list_1_2_cluster$DDE.male)
  gene_list_1_2_cluster$rank_female <- rank(-gene_list_1_2_cluster$DDE.female)

  (gene_list_1_2_cluster$DDE.male > 0)  %>% sum -> x_line
  (gene_list_1_2_cluster$DDE.female > 0)  %>% sum -> y_line

  #(ggplot(gene_list_1_2_cluster, aes(x = rank_female, y = rank_male)) + geom_point() + theme_classic(base_size=20) + ggtitle(clusterID[i])) %>% 
  #		ggrastr::rasterize() %>% print
  #(ggplot(gene_list_1_2_cluster, aes(x = rank_female, y = rank_male)) + geom_hex(binwidth = RRHO2:::defaultStepSize(gene_list1_cluster, gene_list2_cluster)) +
  #              scale_fill_continuous(type = "viridis") + theme_classic(base_size=20) +
  #              geom_vline(xintercept = y_line) + geom_hline(yintercept = x_line) + ggtitle(clusterID[i])) %>%  ggrastr::rasterize() %>% print
  (ggplot(gene_list_1_2_cluster, aes(x = rank_female, y = rank_male)) + geom_density_2d_filled(show.legend = FALSE) + theme_classic(base_size=20) +
                geom_vline(xintercept = y_line) + geom_hline(yintercept = x_line) + ggtitle(clusterID[i])) %>%  ggrastr::rasterize() %>% print
}

dev.off()

names(corr_results) <- c("cluster", "estimate", "statistic", "p.value", "alternative", "method")

write.csv(corr_results, file = "/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_outputs/2_spearman_cluster.csv")

print(summary(warnings()))
sessionInfo()
