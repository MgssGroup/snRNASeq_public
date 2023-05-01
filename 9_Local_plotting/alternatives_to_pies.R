library(tidyverse)
library(randomcoloR)
library(ggvenn)
#library(ggrepel)

base_path <- "C:/Users/Home/OneDrive - McGill University/DifferentialGeneExpression_CombinedData_Interpretation/FilteredDEGS_20220316/"

files <- list(male_broad = c(file_name = "01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct_MaleBroad.csv", 
                          sex = "male", category = "broad"),
           male_subtype = c(file_name = "01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct_MaleSubtype.csv", 
                            sex = "male", category = "subtype"),
           female_broad = c(file_name = "01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct_FemaleBroad.csv",
                            sex = "female", category = "broad"),
           female_subtype = c(file_name = "01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct_FemaleSubtype.csv",
                              sex = "female", category = "subtype"),
           p_broad = c(file_name = "Combination_analysis/9_combination_pvals_all_broad.csv", 
                       sex = "p_combined", category = "broad"),
           p_subtype = c(file_name = "Combination_analysis/9_combination_pvals_all_subtype.csv", 
                         sex = "p_combined", category = "subtype"))

#make color dictionary for consistent colors. 
broads <- c()
clusters <- c()

for(name_file in names(files)) {
  file_to_load <- files[[name_file]]["file_name"]
  filtered_results <- read_csv(file = paste0(base_path, file_to_load)) %>% select(-c(1))
  if(name_file %in% c("p_broad", "p_subtype")) {
    filtered_results %>% filter(signs == 1) %>% mutate(logFC = logFC.male,
                                                  cluster_id = cluster_id.male, .keep = "unused")-> filtered_results
  }
  if(grepl("broad", name_file)) { broads <- union(broads, unique(filtered_results$cluster_id)) }
  if (grepl("subtype", name_file)) { clusters <- union(clusters, unique(filtered_results$cluster_id)) }
}

set.seed(22)
colordict <-list(broad = deframe(bind_cols(broads, distinctColorPalette(length(broads)))), 
                 subtype = deframe(bind_cols(clusters, distinctColorPalette(length(clusters)))))

data <- NULL
pdf(file = paste0(base_path, "DEG_bar_plots.pdf"), onefile = TRUE, height = 14, width = 8)

dataset_genes_list <- list()

for(name_file in names(files)) {
  file_to_load <- files[[name_file]]["file_name"]
  filtered_results <- read_csv(file = paste0(base_path, file_to_load)) %>% select(-c(1))
  if(name_file %in% c("p_broad", "p_subtype")) {
    filtered_results %>% filter(signs == 1) %>% mutate(logFC = logFC.male,
                                                    cluster_id = cluster_id.male,
                                                    p_adj.loc = adjpval,
                                                    .keep = "unused")-> filtered_results
  }
  filtered_results$gene %>% unique %>% head %>% print
  print(name_file)
  dataset_genes_list[[name_file]] <- filtered_results$gene %>% unique
  category <- files[[name_file]]["category"]
  sex <- files[[name_file]]["sex"]
  filtered_results %>% select(cluster_id, logFC, p_adj.loc) -> filtered_results 
  to_add <- setdiff(colordict[[category]] %>% names, filtered_results$cluster_id %>% unique)
  num_to_add <- length(to_add)
  print(to_add)
  filtered_results %>% add_row(cluster_id = to_add, logFC = rep(NA, length(to_add)), p_adj.loc = rep(NA, length(to_add))) %>% 
                      arrange(cluster_id) -> filtered_results
  dp <- ggplot(filtered_results, aes(y = cluster_id, x = logFC, color = p_adj.loc, size = 4)) + geom_point() +
                 geom_vline(xintercept = 0) +
                 guides(size = "none", color = guide_legend(title = "p-value", reverse = TRUE)) +  
                 scale_x_discrete(drop = FALSE) +   
                 theme_classic(base_size = 22) +
                 theme(axis.text.x = element_text(angle = 90))
  print(dp)
  write.csv(dp$data, file = paste0(base_path, name_file, "_DEG_dotplot_source_data.csv" ))
  for(direction in c("up", "down", "all")) {
    data_this <- switch(direction,
                   "up" = filtered_results %>% filter(logFC > 0) %>% group_by(cluster_id) %>% summarise(n()),
                   "down" = , filtered_results %>% filter(logFC < 0) %>% group_by(cluster_id) %>% summarise(n()),
                   "all" = filtered_results %>% select(cluster_id) %>% group_by(cluster_id) %>% summarise(n()))
    data_this %>% mutate(category <- rep(category, n()), sex = rep(sex, n()), direction = rep(direction, n())) -> data_this
    data <- bind_rows(data, data_this)
  }  
}

colnames(data) <- c("cluster_id", "Freq", "category", "sex", "direction")

bar_plotter <- function(data, this_sex, this_category) {
  data %>% filter(sex == this_sex, category == this_category, direction != "all") -> data
  bar_plot <- ggplot(data, aes(y = cluster_id, x=Freq, fill = cluster_id)) +
    geom_bar(stat = "identity") +
    facet_wrap(~direction, ncol = 1) +
    theme_classic(base_size = 22) +
    scale_fill_manual(values = colordict[[this_category]]) +
    ggtitle(paste0(this_sex, "|", this_category)) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90), strip.background = element_rect(linetype = "blank"))
  print(bar_plot)
}


  bar_plotter(data, "male", "broad")
  bar_plotter(data, "male", "subtype")
  bar_plotter(data, "female", "broad")
  bar_plotter(data, "female", "subtype")
  bar_plotter(data, "p_combined", "broad")
  bar_plotter(data, "p_combined", "subtype")
  ggvenn(dataset_genes_list[c("male_subtype", "male_broad", "female_subtype", "female_broad")], text_size = 5) -> p1 
  p1 %>% print
  sink(file = paste0(base_path, "ggvenn_source_data_DEG_overlaps_sex_type.csv"))
  print(dataset_genes_list[c("male_subtype", "male_broad", "female_subtype", "female_broad")]) 
  sink()            
  ggvenn(dataset_genes_list[c("male_subtype", "female_subtype", "p_subtype")], text_size = 5) -> p1
  p1 %>% print
  sink(file = paste0(base_path, "ggvenn_source_data_DEG_overlaps_pcomb_cluster.csv"))
  print(dataset_genes_list[c("male_subtype", "female_subtype", "p_subtype")]) 
  sink()          
  ggvenn(dataset_genes_list[c("male_broad", "female_broad", "p_broad")], text_size = 5) -> p1 
  p1 %>% print
  sink(file = paste0(base_path, "ggvenn_source_data_DEG_overlaps_pcomb_broad.csv"))
  print(dataset_genes_list[c("male_broad", "female_broad", "p_broad")])
  sink()
  
dev.off()

print(colordict)

# $broad
# Ast       Oli       OPC       InN       Mix       Mic 
# "#8AE27B" "#A7D8CE" "#D9CC72" "#BA58D7" "#DF7E84" "#B0A5D7" 
# 
# $subtype
# ExN12_L56    ExN11_L56     ExN3_L46     ExN9_L23    ExN10_L46         Ast1         Ast2     InN3_VIP     InN4_VIP     InN2_SST     InN7_Mix         Oli3 
# "#81A95B"    "#779388"    "#BDE379"    "#57E294"    "#E9764E"    "#E282C2"    "#676893"    "#9DE1A5"    "#645BD3"    "#B6ED4F"    "#C04F72"    "#E5BCED" 
# Oli1         OPC1         OPC2         ExN7     ExN4_L35         End1      InN9_PV      InN1_PV    ExN16_L56        ExN14     ExN2_L23 InN10_ADARB2 
# "#67E65F"    "#DCE9A4"    "#7597E2"    "#72B8E2"    "#DADA3D"    "#EF4C62"    "#DC9A7F"    "#5EE3CB"    "#C090E6"    "#ABECDC"    "#8F38E6"    "#CC67DA" 
# InN8_ADARB2    ExN15_L56     InN5_SST   InN6_LAMP5         Oli2         ExN5          Mix         Mic1    ExN13_L56     ExN1_L24         ExN6 
# "#DAE3C4"    "#CDDBE6"    "#D0BB84"    "#AAA9CD"    "#DD3EE1"    "#E0BEBD"    "#D894AF"    "#ECDF73"    "#E444A8"    "#E1A74D"    "#72D7E2"

print(summary(warnings()))
sessionInfo()

# Summary of (a total of 5) warning messages:
#   1x : Removed 3 rows containing missing values (`geom_point()`).
# 1x : Removed 17 rows containing missing values (`geom_point()`).
# 1x : Removed 1 rows containing missing values (`geom_point()`).
# 1x : Removed 8 rows containing missing values (`geom_point()`).
# 1x : Removed 5 rows containing missing values (`geom_point()`).

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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggvenn_0.1.9        randomcoloR_1.1.0.1 forcats_0.5.2       stringr_1.5.0       dplyr_1.0.10        purrr_1.0.0         readr_2.1.3        
# [8] tidyr_1.2.1         tibble_3.1.8        ggplot2_3.4.0       tidyverse_1.3.2    
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.0    haven_2.5.1         gargle_1.2.1        V8_4.3.0            colorspace_2.0-3    vctrs_0.5.1         generics_0.1.3     
# [8] utf8_1.2.2          rlang_1.0.6         pillar_1.8.1        glue_1.6.2          withr_2.5.0         DBI_1.1.3           bit64_4.0.5        
# [15] dbplyr_2.2.1        modelr_0.1.10       readxl_1.4.1        lifecycle_1.0.3     munsell_0.5.0       gtable_0.3.1        cellranger_1.1.0   
# [22] rvest_1.0.3         labeling_0.4.2      tzdb_0.3.0          parallel_4.2.2      curl_4.3.3          fansi_1.0.3         broom_1.0.2        
# [29] Rcpp_1.0.9          scales_1.2.1        backports_1.4.1     googlesheets4_1.0.1 vroom_1.6.0         jsonlite_1.8.4      farver_2.1.1       
# [36] bit_4.0.5           fs_1.5.2            hms_1.1.2           stringi_1.7.8       Rtsne_0.16          cli_3.5.0           tools_4.2.2        
# [43] magrittr_2.0.3      cluster_2.1.4       crayon_1.5.2        pkgconfig_2.0.3     ellipsis_0.3.2      xml2_1.3.3          reprex_2.0.2       
# [50] googledrive_2.0.0   lubridate_1.9.0     timechange_0.1.1    assertthat_0.2.1    httr_1.4.4          rstudioapi_0.14     R6_2.5.1           
# [57] compiler_4.2.2 