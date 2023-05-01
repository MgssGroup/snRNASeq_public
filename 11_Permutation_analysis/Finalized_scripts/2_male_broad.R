library(tidyverse)
library(limma)
library(muscat)
library(Seurat)
library(SingleCellExperiment)
library(BiocParallel)

print("Set up parallel processing for muscat based on https://github.com/HelenaLC/muscat/issues/70")
registered() # I have MulticoreParam capabilities
options(MulticoreParam=quote(MulticoreParam()))
param <- MulticoreParam(log = FALSE)

print("Set input and output directory")
input_dir <-
  "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/"
output_dir <-
  "/home/malosree/scratch/Permutation_outputs/2_male_broad/"

print("Read data")
harmonized_object <- readRDS(file = paste0(input_dir, "1_enhanced_harmonized_object.Rds"))

print("Subset data")
MC_raw_clustered_Male <-
  subset(harmonized_object, subset = (Sample != c("M24_2") & Sex == "Male"))
rm(harmonized_object)

original_case_control <- unique.data.frame(MC_raw_clustered_Male@meta.data[c("Sample",
                           "Condition", "Batch")]) %>% arrange(Batch, Sample)

args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args)

permute_labels <- function(i, x, y) {
        set.seed(i)
        x$Condition <- sample(x$Condition, size = length(x$Condition), replace = FALSE)
        x$Batch <- y$Batch
        x
}

#for(i in 1:3) {

	print(paste("Add permuted case-control labels", i))
        permuted_case_control <- original_case_control
	
	permuted_case_control %>% group_by(Batch) %>% group_map(~permute_labels(i, .x, .y)) %>% bind_rows() -> permuted_case_control

        #print(permuted_case_control)
        for (sample in permuted_case_control$Sample) {
		MC_raw_clustered_Male@meta.data$Condition[MC_raw_clustered_Male@meta.data$Sample == sample] <-	
			rep(permuted_case_control$Condition[permuted_case_control$Sample == sample], 
				sum(MC_raw_clustered_Male@meta.data$Sample == sample))
	}
        #print(unique.data.frame(MC_raw_clustered_Male@meta.data[c("Sample",
        #                   "Condition")]))
	
	print(paste("Make and preprocess SingleCellExperiment object", i))
	sce <- as.SingleCellExperiment(MC_raw_clustered_Male)

	rm(MC_raw_clustered_Male)

	sce <- sce[rowSums(counts(sce) > 0) > 0, ]

	sce <- prepSCE(
		sce,
		kid = "Broad",
		sid = "Sample",
		gid = "Condition",
		drop = FALSE)
  
	print(paste("Make design matrix", i))
	#contrast####
	#design & contrast matrix
	temp <- colData(sce)
	rownames(temp) <- NULL
	anno <-
		unique.data.frame(temp[c("sample_id",
			"group_id",
			"Sex",
			"Batch",
			"Age",
			"pH",
			"PMI")])

	design <- model.matrix(~0+group_id+Batch+Age+PMI+pH, data = anno)

	dimnames(design) <-
		list(anno$sample_id, levels = make.names(colnames(design)))
	print(design)

	print(paste("Make contrast", i))
	# Making a contrast for each pairwise comparison.
	contrast <-
	makeContrasts(
		Case_vs_Control = (group_idCase - group_idControl),
			levels = make.names(colnames(design)))
  
	print(paste("Aggregate data", i))
	#aggData####
	pb <- aggregateData(
		sce,
		assay = "counts",
		fun = "sum",
		by = c("cluster_id", "sample_id"),
		BPPARAM = param)
  
	print(paste("Perform pseudobulk differential expression", i))
	#pbDS####
	# run DS analysis c("edgeR", "DESeq2", "limma-trend", "limma-voom")
	res <- pbDS(
		pb,
		design = design,
		min_cells = 10,
		contrast = contrast,
		filter = "none",
		method = "edgeR", 
		BPPARAM = param)
  
	print(paste("Format and save results", i))  
	#resDS####
	#output_master# append CPMs & expression frequencies
	try(res <-
		resDS(sce,
		res,
		cpm = TRUE,
		digits = 10,
		bind = "row"))
 
	# tidy format
        #Method_name####
        write_csv(res, file = paste0(output_dir,i,
                 "_unfiltered_male_broad.csv"))
 
#}

print(summary(warnings()))
sessionInfo()
