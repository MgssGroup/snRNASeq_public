library(tidyverse)
library(ggrepel)

base_path <- "C:/Users/Home/OneDrive - McGill University/DifferentialGeneExpression_CombinedData_Interpretation/FilteredDEGS_20220316/Combination_analysis/"

files <- c(p_broad = "9_combination_pvals_all_broad.csv",
           p_subtype = "9_combination_pvals_all_subtype.csv",
           model_broad = "Combined_model_filteredBroad.csv",
           model_subype = "Combined_model_filteredCluster.csv")

print("Make pie plots of numbers of DEGs per cluster after filtering")
dirs <- c("Down", "Up")
pie_plots <- list()

for(file_name in names(files)) {
  file_to_load <- files[[file_name]]
	filtered_results <- read_csv(file = paste0(base_path, file_to_load)) %>% select(-c(1)) 
	if(file_name %in% c("p_broad", "p_subtype")) filtered_results %>% filter(signs == 1) %>% mutate(logFC = logFC.male, 
	                                                cluster_id = cluster_id.male, .keep = "unused")-> filtered_results
	for(direction in c("up", "down", "all")) {
	  print(direction)
	  data <- switch(direction,
	                 "up" = filtered_results %>% filter(logFC > 0) %>% select(cluster_id) %>% unlist() %>% table() %>% as.data.frame(),
	                 "down" = , filtered_results %>% filter(logFC < 0) %>% select(cluster_id) %>% unlist() %>% table() %>% as.data.frame(),
	                 "all" = filtered_results %>% select(cluster_id) %>% unlist() %>% table() %>% as.data.frame()
	                 #"combined" = filtered_results %>% mutate(dir = map_chr(logFC, ~dirs[sign(.x)/2+1.5])) %>%
	                 #  mutate(cluster_id = paste0(dir, "_", cluster_id)) %>%
	                 #  select(cluster_id) %>% unlist() %>% table() %>% as.data.frame()
	                 )
	  print(dim(data))
	  colnames(data) <- c("cluster_id", "Freq")
	  data <- data %>%
	    mutate(csum = rev(cumsum(rev(Freq))),
	           pos = Freq/2 + lead(csum, 1),
	           pos = if_else(is.na(pos), Freq/2, pos)) %>%
	    mutate(cluster_id = paste(cluster_id, Freq, sep = " - "))
	  pie_plots[[paste0(direction, "|" , file_name)]] <- ggplot(data, aes(x = "", y=Freq, fill=cluster_id)) +
	    geom_bar(stat="identity", width=1) +
	    coord_polar("y", start=0) +
	    theme_void(base_size = 22) +
	    #geom_label_repel(aes(y = pos), color = "black", size=6, show.legend = FALSE, nudge_x = 1) +
	    ggtitle(paste0(direction, "|" , file_name))
	}
}

p_comb <- read_csv(file = "C:/Users/Home/OneDrive - McGill University/DifferentialGeneExpression_CombinedData_Interpretation/FilteredDEGS_20220316/Combination_analysis/9_combination_pvals_all_subtype.csv") %>% filter(signs == 1)
female <- read_csv(file = "C:/Users/Home/OneDrive - McGill University/DifferentialGeneExpression_CombinedData_Interpretation/FilteredDEGS_20220316/01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct_FemaleSubtype.csv")
male <- read_csv(file = "C:/Users/Home/OneDrive - McGill University/DifferentialGeneExpression_CombinedData_Interpretation/FilteredDEGS_20220316/01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct_MaleSubtype.csv")
table(p_comb$cluster_id.female) %>% as.data.frame() -> p_comb
table(female$cluster_id) %>% as.data.frame() -> female
table(male$cluster_id) %>% as.data.frame() -> male

bind_rows(list(p_comb = p_comb, male = male, female = female), .id = "Analysis") -> plotted

pdf(file = paste0(base_path, "combined_pie_charts.pdf"), onefile = TRUE, height = 7, width = 10)
        lapply(pie_plots, print)
      plt <- (ggplot(data = plotted, aes(x = Var1, y = Freq, color = Analysis)) + geom_point(size = 3)+ 
        theme_classic(base_size = 20) + 
        theme(axis.text.x = element_text(angle = 45, size = 16, hjust = 1), axis.title.x = element_blank(), 
              axis.title.y = element_blank(), legend.position = "bottom"))
      print(plt)
      write.csv(plt$data, file = paste0(base_path, "combined_DEGS_scatter_source_data.csv"))
dev.off()

p_comb_broad <- read_csv(file = "C:/Users/Home/OneDrive - McGill University/DifferentialGeneExpression_CombinedData_Interpretation/FilteredDEGS_20220316/Combination_analysis/9_combination_pvals_all_broad.csv") %>% filter(signs == 1)
female_broad <- read_csv(file = "C:/Users/Home/OneDrive - McGill University/DifferentialGeneExpression_CombinedData_Interpretation/FilteredDEGS_20220316/01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct_FemaleBroad.csv")
male_broad <- read_csv(file = "C:/Users/Home//OneDrive - McGill University/DifferentialGeneExpression_CombinedData_Interpretation/FilteredDEGS_20220316/01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct_MaleBroad.csv")

clusters <- reduce(c(unique(p_comb_broad$cluster_id.female), unique(male_broad$cluster_id), unique(female_broad$cluster_id)), union)

num_overlap <- NULL

for(cluster in clusters) {
  common <- length(intersect(p_comb_broad$gene[p_comb_broad$cluster_id.female == cluster], male_broad$gene[male_broad$cluster_id == cluster]))
  num_overlap <- rbind(num_overlap, c(cluster, "male", common, length(p_comb_broad$gene[p_comb_broad$cluster_id.female == cluster]), length(male_broad$gene[male_broad$cluster_id == cluster])))
  common <- length(intersect(p_comb_broad$gene[p_comb_broad$cluster_id.female == cluster], female_broad$gene[female_broad$cluster_id == cluster]))
  num_overlap <- rbind(num_overlap, c(cluster, "female", common, length(p_comb_broad$gene[p_comb_broad$cluster_id.female == cluster]), length(female_broad$gene[female_broad$cluster_id == cluster])))
}

colnames(num_overlap) <- c("Cluster", "Sex", "Common", "P-comb", "Sex-specific")

write.csv(num_overlap, file = paste0(base_path, "broad_cluster_by_cluster_overlap.csv"))

print(summary(warnings()))
# Length  Class   Mode 
# 0   NULL   NULL 
sessionInfo()
# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_Canada.utf8  LC_CTYPE=English_Canada.utf8    LC_MONETARY=English_Canada.utf8 LC_NUMERIC=C                   
# [5] LC_TIME=English_Canada.utf8    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggrepel_0.9.3   forcats_0.5.2   stringr_1.5.0   dplyr_1.0.10    purrr_1.0.0     readr_2.1.3     tidyr_1.2.1    
# [8] tibble_3.1.8    ggplot2_3.4.0   tidyverse_1.3.2
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.0    haven_2.5.1         gargle_1.2.1        colorspace_2.0-3    vctrs_0.5.1         generics_0.1.3     
# [7] utf8_1.2.2          rlang_1.0.6         pillar_1.8.1        glue_1.6.2          withr_2.5.0         DBI_1.1.3          
# [13] bit64_4.0.5         dbplyr_2.2.1        modelr_0.1.10       readxl_1.4.1        lifecycle_1.0.3     munsell_0.5.0      
# [19] gtable_0.3.1        cellranger_1.1.0    rvest_1.0.3         labeling_0.4.2      tzdb_0.3.0          parallel_4.2.2     
# [25] fansi_1.0.3         broom_1.0.2         Rcpp_1.0.9          scales_1.2.1        backports_1.4.1     googlesheets4_1.0.1
# [31] vroom_1.6.0         jsonlite_1.8.4      farver_2.1.1        fs_1.5.2            bit_4.0.5           hms_1.1.2          
# [37] stringi_1.7.8       grid_4.2.2          cli_3.5.0           tools_4.2.2         magrittr_2.0.3      crayon_1.5.2       
# [43] pkgconfig_2.0.3     ellipsis_0.3.2      xml2_1.3.3          reprex_2.0.2        googledrive_2.0.0   lubridate_1.9.0    
# [49] timechange_0.1.1    assertthat_0.2.1    httr_1.4.4          rstudioapi_0.14     R6_2.5.1            compiler_4.2.2