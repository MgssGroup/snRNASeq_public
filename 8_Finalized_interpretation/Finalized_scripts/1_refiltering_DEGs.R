library(tidyverse)
library(scater)

base_path <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/"

directories <- c("male_broad" = "04_Male_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Male_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/", 
		"male_subtype" = "Mar9_2022_updated_res/", 
		"female_broad" = "03_Female_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/", 
		"female_subtype" = "03_Female_dataset_analysis/01_pseudobulk/02_cell_subtype/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_subtype/")

metadata <- read_csv(file = paste0(base_path, "Male_female_metadata_combined.csv"))

case_subjects <- metadata %>% filter(Condition == "Case") %>% select(Sample) %>% unlist() %>% paste0(".cpm")
print(case_subjects)
control_subjects <- metadata %>% filter(Condition == "Control") %>% select(Sample) %>% unlist() %>% paste0(".cpm")
print(control_subjects)

#Log version added 2022.04.18
count_outliers <- function(x, log = FALSE) {isOutlier(x[!is.na(x)], nmads =5, log = log) %>% sum }

for(directory_name in names(directories)) {
	print("Read in data")
	
	print(directory_name)

	directory <- directories[[directory_name]]	
	if(directory_name == "male_broad") {
                results <- readRDS(file = paste0(base_path, directory, "res.rds"))
        } else {
                res <- readRDS(paste0(base_path, directory, "res.rds"))
                sce <- readRDS(paste0(base_path, directory, "sce.rds"))
                results <- muscat::resDS(sce, res, bind = "row", cpm = TRUE, frq = FALSE, digits = 10)
        }

	print("Calculate additional parameters to add to the DEG table")
	results %>% mutate(NumNonZero = rowSums(across(ends_with(".cpm"), ~(.x>0)), na.rm = TRUE), 
			NumCaseExcluded = rowSums(across(any_of(case_subjects), is.na)), 
			NumControlExcluded = rowSums(across(any_of(control_subjects), is.na))) %>%
			mutate(Greater3 = NumNonZero > 2) -> results

	print("Filter based on p_value, FC, and number of subjects with gene detected")
	results %>% filter(p_adj.loc < 0.05, abs(logFC) > log2(1.1), Greater3 == TRUE) -> filtered_results

	filtered_results %>% rowwise() %>% mutate(NumOutliers = count_outliers(across(ends_with(".cpm"))), 
						NumOutliersLog = count_outliers(across(ends_with(".cpm")), log = TRUE)) -> filtered_results

	print("Write results")
	write.csv(results, file = paste0(base_path,
                                directory,
                                "01_Table_raw_masterfile_cpm_Freq_DS_stats_greater3.csv"))
	write.csv(filtered_results, file = paste0(base_path,
                                directory,
                                "01_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct.csv"))
}

print(summary(warnings()))
sessionInfo()

