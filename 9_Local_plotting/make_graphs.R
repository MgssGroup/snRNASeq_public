library(tidyverse)
library(igraph)
library(ggvenn)

setwd("C:/Users/Home/OneDrive - McGill University/DifferentialGeneExpression_CombinedData_Interpretation/FGSEA_download/Collapsed_pathways")

to_remove_InN1_PV <- c("REACTOME_ELEVATION_OF_CYTOSOLIC_CA2_LEVELS",
                       "REACTOME_OLFACTORY_SIGNALING_PATHWAY", 
                       "REACTOME_MECP2_REGULATES_NEURONAL_RECEPTORS_AND_CHANNELS")

clusters <- c("Mic1" =  "7_collapsed_pathways_female_subtype_Mic1_CP_REACTOME",
              "InN1" = "7_collapsed_pathways_female_subtype_InN1_PV_CP_REACTOME",
              "InN9" = "7_collapsed_pathways_female_subtype_InN9_PV_CP_REACTOME")

clean_pathway <- function(pathway, type) {
  str_replace(pathway, pattern = fixed(type), replacement = "") %>% str_replace_all(pattern = "_", replacement = " ") %>% str_to_title() -> pathway
}

for(cluster in names(clusters)) {
  df <- read_csv(file = paste0(clusters[[cluster]], ".csv")) %>% setNames(nm = c("Pathway", "Parent"))
  df %>% mutate(Pathway = map_chr(Pathway, ~clean_pathway(.x, "REACTOME_" ))) %>% mutate(Parent = map_chr(Parent, ~clean_pathway(.x, "REACTOME_"))) -> df
  df %>% drop_na %>% as.data.frame %>% as.matrix %>% graph_from_edgelist -> df_graph
  df_graph %>% print
  write_graph(df_graph, file = paste0(clusters[[cluster]], ".gml"), format = "gml")
}

df <- list()

df[["MEturquoise"]] <- read_csv("GOreactome_turquoise.csv") %>% select(ID) %>% unlist
df[["InN_PV_down"]] <- read_csv("7_collapsed_pathways_female_subtype_InN1_PV_CP_REACTOME_reactome_ids.csv") %>% 
        filter(!(gs_name %in% to_remove_InN1_PV)) %>%
        select(gs_exact_source) %>% unlist
df[["InN_PV_down"]] <- union(df[["InN_PV_down"]], read_csv("7_collapsed_pathways_female_subtype_InN9_PV_CP_REACTOME_reactome_ids.csv") %>% 
 select(gs_exact_source) %>% unlist)

sink(file = "InPV_reactomeoverlap_soruce_data.csv")
print(df)
sink()

pdf(file = "InPV_reactomeoverlap.pdf")
ggvenn(data = df, text_size = 5)
dev.off()

print(summary(warnings()))
# 1 identical warnings:
#   One or more parsing issues, call `problems()` on your data frame for details, e.g.:
#   dat <- vroom(...)
# problems(dat)

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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] igraph_1.3.5        ggvenn_0.1.9        randomcoloR_1.1.0.1 forcats_0.5.2       stringr_1.5.0       dplyr_1.0.10       
# [7] purrr_1.0.0         readr_2.1.3         tidyr_1.2.1         tibble_3.1.8        ggplot2_3.4.0       tidyverse_1.3.2    
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.0    haven_2.5.1         gargle_1.2.1        V8_4.3.0            colorspace_2.0-3    vctrs_0.5.1        
# [7] generics_0.1.3      utf8_1.2.2          rlang_1.0.6         pillar_1.8.1        glue_1.6.2          withr_2.5.0        
# [13] DBI_1.1.3           bit64_4.0.5         dbplyr_2.2.1        modelr_0.1.10       readxl_1.4.1        lifecycle_1.0.3    
# [19] munsell_0.5.0       gtable_0.3.1        cellranger_1.1.0    rvest_1.0.3         labeling_0.4.2      tzdb_0.3.0         
# [25] parallel_4.2.2      curl_4.3.3          fansi_1.0.3         broom_1.0.2         Rcpp_1.0.9          scales_1.2.1       
# [31] backports_1.4.1     googlesheets4_1.0.1 vroom_1.6.0         jsonlite_1.8.4      farver_2.1.1        bit_4.0.5          
# [37] fs_1.5.2            hms_1.1.2           stringi_1.7.8       Rtsne_0.16          cli_3.5.0           tools_4.2.2        
# [43] magrittr_2.0.3      cluster_2.1.4       crayon_1.5.2        pkgconfig_2.0.3     ellipsis_0.3.2      xml2_1.3.3         
# [49] reprex_2.0.2        googledrive_2.0.0   lubridate_1.9.0     timechange_0.1.1    assertthat_0.2.1    httr_1.4.4         
# [55] rstudioapi_0.14     R6_2.5.1            compiler_4.2.2     