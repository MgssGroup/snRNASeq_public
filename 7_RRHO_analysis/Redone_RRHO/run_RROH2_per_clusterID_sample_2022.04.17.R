closeAllConnections()
rm(list = ls())
#library(GGally)
library(ggpubr)
library(patchwork)
library(RRHO2)
library(tidyverse)
library(cowplot)
#library(ggeasy)
library(gplots)

summary_table <- data.frame()

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
clusterID <- intersect(unique(res_female_raw$cluster_id), unique(res_male_raw$cluster_id)) %>% sort()
for(method_type in c("hyper")) {
#c("fisher", "hyper")) {
for(correction in c("none", "BY")) {
#c("none", "BY", "BH")) {
for(scaling in c(TRUE)) {
#c(TRUE, FALSE)) {
pdf(paste0("Redone_RRHO/Figure_",method_type, "_", correction, "_", scaling, "_RRHO2_correlation_cell_Broad.pdf"),
    width = 6,
    height = 4)
for (i in 1:length(clusterID)) {
  #i = 2
  print(paste0("Running clusterID = ", clusterID[i]))
  #i =1
  gene_list1_cluster <- gene_list1[grep(paste0("_", clusterID[i]),gene_list1$Genes),]
  gene_list2_cluster <- gene_list2[grep(paste0("_", clusterID[i]),gene_list2$Genes),]
  gene_list2_cluster <- gene_list2_cluster[match(gene_list1_cluster$Genes, gene_list2_cluster$Genes), ]

  RRHO_obj<- NULL
  try(RRHO_obj <-
        RRHO2_initialize(gene_list2_cluster,
                         gene_list1_cluster,
                         labels = c(paste0("Female_",clusterID[i]),
				#"_Batch_Age_PMI_pH"), 
				paste0("Male_",clusterID[i])),
				#"_Batch_Age_PMI_pH")),
                         log10.ind = scaling, method = method_type, multipleTesting = correction))
  
  print("Genes contributing to overlaps")
  head(RRHO_obj$genelist_dd$gene_list_overlap_dd) %>% print()
  head(RRHO_obj$genelist_uu$gene_list_overlap_uu) %>% print()
  head(RRHO_obj$genelist_du$gene_list_overlap_du) %>% print()
  head(RRHO_obj$genelist_ud$gene_list_overlap_ud) %>% print()
  max(RRHO_obj$hypermat, na.rm = TRUE) %>% print()	

  write.csv(RRHO_obj$hypermat, file = paste0(
	"Redone_RRHO/Source_data_RRHO/Source_data", clusterID[i],
	method_type, "_", correction, "_", scaling, "_RRHO2_correlation_cell_Broad.csv")) 
  
  summary_table <- rbind(summary_table, c(clusterID[i], method_type, correction, scaling, max(RRHO_obj$hypermat, na.rm = TRUE)))	
	
  try(RRHO2_heatmap(RRHO_obj))
  if(method_type == "hyper" & scaling == TRUE & correction == "none")  {
	try(RRHO2_heatmap(RRHO_obj, maximum = 50))
	#try(RRHO2_heatmap(RRHO_obj, maximum = 100))
  }
  try(RRHO2_vennDiagram(RRHO_obj, "dd"))
  try(RRHO2_vennDiagram(RRHO_obj, "uu"))
  try(RRHO2_vennDiagram(RRHO_obj, "du"))
  try(RRHO2_vennDiagram(RRHO_obj, "ud"))

}
dev.off()
}  
}		
}


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

clusterID <- intersect(unique(res_female_raw$cluster_id), unique(res_male_raw$cluster_id)) %>% sort()
for(method_type in c("hyper")) {
#c("fisher", "hyper")) {
for(correction in c("none", "BY")) {
#c("none", "BY", "BH")) {
for(scaling in c(TRUE)) {
#c(TRUE, FALSE)) {
  pdf(paste0("Redone_RRHO/Figure_",method_type, "_", correction, "_", scaling, "_RRHO2_correlation_cell_subtype.pdf"),
    width = 6,
    height = 4)
for (i in 1:length(clusterID)) {
  #i = 2
  print(paste0("Running clusterID = ", clusterID[i]))
  #i =1
  gene_list1_cluster <- gene_list1[grep(paste0("_", clusterID[i]),gene_list1$Genes),]
  gene_list2_cluster <- gene_list2[grep(paste0("_", clusterID[i]),gene_list2$Genes),]
  gene_list2_cluster <- gene_list2_cluster[match(gene_list1_cluster$Genes, gene_list2_cluster$Genes), ]

  try(RRHO_obj <-
        RRHO2_initialize(gene_list2_cluster,
                         gene_list1_cluster,
                         labels = c(paste0("Female_",clusterID[i]
				#,"_Batch_Age_PMI_pH"
				), paste0("Male_",clusterID[i]
				#,"_Batch_Age_PMI_pH"
				)),
                         log10.ind = scaling, method = method_type, multipleTesting = correction))
  
  print("Genes contributing to overlaps")
  head(RRHO_obj$genelist_dd$gene_list_overlap_dd) %>% print()
  head(RRHO_obj$genelist_uu$gene_list_overlap_uu) %>% print()
  head(RRHO_obj$genelist_du$gene_list_overlap_du) %>% print()
  head(RRHO_obj$genelist_ud$gene_list_overlap_ud) %>% print()
  max(RRHO_obj$hypermat, na.rm = TRUE) %>% print

  write.csv(RRHO_obj$hypermat, file = paste0(
	"Redone_RRHO/Source_data_RRHO/Source_data_", clusterID[i],
	method_type, "_", correction, "_", scaling, "_RRHO2_correlation_cell_subtype.csv"))

  summary_table <- rbind(summary_table, c(clusterID[i], method_type, correction, scaling, max(RRHO_obj$hypermat, na.rm = TRUE)))

  try(RRHO2_heatmap(RRHO_obj))
  if(method_type == "hyper" & scaling == TRUE & correction == "none")  {  
	try(RRHO2_heatmap(RRHO_obj, maximum = 25))
	#try(RRHO2_heatmap(RRHO_obj, maximum = 100))	
  }
  try(RRHO2_vennDiagram(RRHO_obj, "dd"))
  try(RRHO2_vennDiagram(RRHO_obj, "uu"))
  try(RRHO2_vennDiagram(RRHO_obj, "du"))
  try(RRHO2_vennDiagram(RRHO_obj, "ud"))

}
dev.off()
}
}
}

colnames(summary_table) <- c("Cluster", "Method", "Correction", "Scaling", "Max_p")

write.csv(summary_table, "Redone_RRHO/summary_table_RRHO.csv")

print(summary(warnings()))
sessionInfo()
