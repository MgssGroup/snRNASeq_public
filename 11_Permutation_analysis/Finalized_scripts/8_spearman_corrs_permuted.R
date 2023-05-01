library(tidyverse)

all_results_broad <- tibble()

all_results_cluster <- tibble()

#Broad#####

setwd("/home/malosree/projects/def-cnagy/malosree/Permutation_outputs/")

for(iter in 1:100) { 
  res_male_raw <-
    read_csv(file = paste0("2_male_broad/", iter, "_unfiltered_male_broad.csv")) %>% 
    as.data.frame()
  res_male_raw$ID <-
    paste0(res_male_raw$gene, "_", res_male_raw$cluster_id)
  res_male <- select(res_male_raw, ID, p_val, logFC)	
  res_male$Score  <-  -log10(res_male$p_val)
  res_male$Score <- res_male$Score * (res_male$logFC)
  
  gene_list1 <-
    data.frame(Genes = res_male$ID, DDE = res_male$Score)
  
  res_female_raw <-
    read_csv(paste0("4_female_broad/", iter, "_unfiltered_female_broad.csv")) %>% 
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
  
  clusterID <- intersect(unique(res_female_raw$cluster_id), unique(res_male_raw$cluster_id)) %>% sort()
  for (i in 1:length(clusterID)) {
    #i = 2
    print(paste0("Running clusterID = ", clusterID[i], ", iteration: ", iter))
    #i =1
    gene_list1_cluster <- gene_list1[grep(paste0("_", clusterID[i]),gene_list1$Genes),]
    gene_list2_cluster <- gene_list2[grep(paste0("_", clusterID[i]),gene_list2$Genes),]
    gene_list2_cluster <- gene_list2_cluster[match(gene_list1_cluster$Genes, gene_list2_cluster$Genes), ]
    
    gene_list_1_2_cluster <- inner_join(gene_list1_cluster, gene_list2_cluster, by = c("Genes"), suffix = c(".male", ".female"))
    
    cor_res <- cor.test(gene_list_1_2_cluster$DDE.female, gene_list_1_2_cluster$DDE.male, method = "spearman")	
    cor_res <- c(clusterID[i], iter, cor_res[c("estimate", "statistic", "p.value", "alternative", "method")])	
    names(cor_res) <- c("cluster", "iteration", "estimate", "statistic", "p.value", "alternative", "method")
    
    corr_results <- rbind(corr_results, cor_res)	
  
  }
  
  names(corr_results) <- c("cluster", "iteration", "estimate", "statistic", "p.value", "alternative", "method")
  
  all_results_broad <- bind_rows(all_results_broad, corr_results)
  
  #cell_subtype####
  
  res_male_raw <-
    read_csv(file = paste0("1_male_subtype/1_male_subtype/", iter, "_unfiltered_male_subtype.csv")) %>% as.data.frame()
  res_male_raw$ID <-
    paste0(res_male_raw$gene, "_", res_male_raw$cluster_id)
  res_male <- select(res_male_raw, ID, p_val, logFC)
  res_male$Score  <- -log10(res_male$p_val)
  res_male$Score <- res_male$Score * (res_male$logFC)
  
  gene_list1 <-
    data.frame(Genes = res_male$ID, DDE = res_male$Score)
  
  res_female_raw <-
    read_csv(file = paste0("3_female_subtype/", iter, "_unfiltered_female_subtype.csv")) %>% as.data.frame()
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
    cor_res <- c(clusterID[i], iter, cor_res[c("estimate", "statistic", "p.value", "alternative", "method")])
    names(cor_res) <- c("cluster", "iteration", "estimate", "statistic", "p.value", "alternative", "method")
    
    corr_results <- rbind(corr_results, cor_res)
    
  }
  
  names(corr_results) <- c("cluster", "iteration", "estimate", "statistic", "p.value", "alternative", "method")
  
  all_results_cluster <- bind_rows(all_results_cluster, corr_results)
  
  
}  

write.csv(all_results_broad, "8_spearman_corrs_permuted_broad.csv")
write.csv(all_results_cluster, "8_spearman_corrs_permuted_cluster.csv")

original_corrs_broad <- read_csv(file = "/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_outputs/2_spearman_broad.csv")
original_corrs_subtype <- read_csv(file = "/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_outputs/2_spearman_cluster.csv")

all_results_cluster$cluster <- recode(all_results_cluster$cluster, "InN8_Mix" = "InN8_ADARB2")

fractions_below <- tibble()

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Finalized_outputs/8_spearman_corr_distributions.pdf", 
	onefile = TRUE)
	for(clus in original_corrs_broad$cluster) {
		par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
		MASS::truehist(all_results_broad$estimate[all_results_broad$cluster == clus], 
			xlab = paste0(clus, " Spearman corr dist"), prob = FALSE)
		abline(v = original_corrs_broad$estimate[original_corrs_broad$cluster == clus])
		fractions_below <- bind_rows(fractions_below, tibble(cluster = clus, 
			fraction = cume_dist(
				c(original_corrs_broad$estimate[original_corrs_broad$cluster == clus], 
					all_results_broad$estimate[all_results_broad$cluster == clus]))[1]))
		write.csv(all_results_broad$estimate[all_results_broad$cluster == clus], file = paste0(
		"/home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Finalized_outputs/8_source_data_", 
		clus, "spearmans_permuted.csv"))
	}
	for(clus in original_corrs_subtype$cluster) {
		print(clus)
                par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
                MASS::truehist(all_results_cluster$estimate[all_results_cluster$cluster == clus],
                        xlab = paste0(clus, " Spearman corr dist"), prob = FALSE)
                abline(v = original_corrs_subtype$estimate[original_corrs_subtype$cluster == clus])
		fractions_below <- bind_rows(fractions_below, tibble(cluster = clus,
                        fraction = cume_dist(
                                c(original_corrs_subtype$estimate[original_corrs_subtype$cluster == clus],
                                        all_results_cluster$estimate[all_results_cluster$cluster == clus]))[1]))        
}
dev.off()

write.csv(fractions_below, file = "/home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Finalized_outputs/8_spearman_fractions_permutations.csv")

print(summary(warnings()))
sessionInfo()
