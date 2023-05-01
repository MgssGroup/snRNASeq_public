library(Seurat)
library(harmony)
library(tidyverse)
library(psych)

harmonized_object <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/3_harmonized_object.Rds")

#Set seed (https://github.com/immunogenomics/harmony/issues/56)
set.seed(23)

#Reharmonize
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/Finalized_outputs/5_seeded_harmony_convergence_plot.pdf", onefile =TRUE)
harmonized_object <- RunHarmony(harmonized_object, reduction.save = "Harmony_Seeded_Batch_Sample_Chemistry", 
		group.by.vars = c("Batch", "Sample", "Chemistry"), kmeans_init_iter_max = 5000,
                plot_convergence = TRUE)
print(summary(warnings()))
dev.off()

runs <- list()

#Check if seeding actually makes things constant.
for(i in 1:3) {
set.seed(23)
temp <- RunHarmony(harmonized_object, reduction.save = "Harmony_Seeded_Batch_Sample_Chemistry", group.by.vars = c("Batch", "Sample", "Chemistry"), kmeans_init_iter_max = 5000,
                plot_convergence = FALSE)
runs[i] <- max(abs(
	harmonized_object@reductions$Harmony_Seeded_Batch_Sample_Chemistry@cell.embeddings - temp@reductions$Harmony_Seeded_Batch_Sample_Chemistry@cell.embeddings))
}
print("Value of runs:")
print(runs)
rm(temp)

#check correspondence to previous run
print("Correlation of each Harmony component to previous run without seed.")

dim_corrs <- vector()
for(i in 1:100) {
	dim_corrs <- rbind(dim_corrs, c(i,
			cor(harmonized_object@reductions$Harmony_Batch_Sample_Chemistry@cell.embeddings[,i], harmonized_object@reductions$Harmony_Seeded_Batch_Sample_Chemistry@cell.embeddings[,i])))
}
write.csv(dim_corrs, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/Finalized_outputs/5_Harmony_seeded_unseeded_corrs.csv")

#just in case
harmonized_object <- UpdateSeuratObject(harmonized_object)

#Elbow plot
print("Plotting elbow")
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/Finalized_outputs/5_Seeded_Harmony_Elbow_plot.pdf", onefile = TRUE)
ElbowPlot(harmonized_object, reduction = "Harmony_Seeded_Batch_Sample_Chemistry", ndims = 100)
dev.off()

#Run UMAP
print("Running UMAP")
harmonized_object <- RunUMAP(harmonized_object, dims = 1:100, reduction= "Harmony_Seeded_Batch_Sample_Chemistry",
                                        reduction.name = "UMAPHarmonySeededBatchSampleChemistry",
                                        reduction.key = "UMAPHarmonySeededBatchSampleChemistry")

#Check UMAP variation by categorial variable numerically
print("UMAP (all 100 dimensions) lms")

variables <- c("Batch", "Chemistry", "Sample")

var_results <- vector()

cell.embeddings <- Embeddings(harmonized_object, reduction = "UMAPHarmonySeededBatchSampleChemistry")
        results_1 <- vector()
        results_2 <- vector()
        for(variable in variables) {
                print(variable)
                results_1 <- c(results_1, summary(lm(cell.embeddings[,1]~harmonized_object@meta.data[, variable]))$r.squared)
                results_2 <- c(results_2, summary(lm(cell.embeddings[,2]~harmonized_object@meta.data[, variable]))$r.squared)
                if(length(results_2) == length(variables)) var_results <- rbind(var_results, results_1, results_2)
        }

rownames(var_results) <- paste("UMAP_Harmony_Seeded_Batch_Sample_Chemistry", 1:2, sep = "_")
colnames(var_results) <- variables

write.csv(var_results, file =
        "/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/Finalized_outputs/5_UMAP_Harmony_Seeded_categorical_vars_lm_all_100.csv")

#Redo some visualizations

#Plot how much UMAP corresponds to batch and other categorical variables.
print("Plot categoricals")

plots_list <- list()

categoricals  <- c("Batch", "Chemistry", "Sex", "Sample")

for(categorical in categoricals) {
       plots_list[[categorical]] <- DimPlot(harmonized_object, reduction= "UMAPHarmonySeededBatchSampleChemistry",
                                                                         group.by = categorical)
}

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/Finalized_outputs/5_UMAPS_Harmony_Seeded_categorical.pdf",
        onefile = TRUE, width = 12)
plots_list
dev.off()

#Plot how much UMAP corresponds to technical continuous variables

#calcuclate ribosomal percentage (file downloaded for group 1054 from HGNC, 2021.06.19)

library(readr)

ribo_proteins <- read_tsv(file = "/home/malosree/scratch/ribo_proteins_hgnc_cgi_1054.tsv",
                                col_names = TRUE, col_types = "ccccccccccccc")
spec(ribo_proteins)
ribo_proteins <- ribo_proteins$`Approved symbol`
ribo_proteins <- intersect(ribo_proteins, rownames(harmonized_object))

harmonized_object <- PercentageFeatureSet(harmonized_object, features = ribo_proteins, col.name = "percent.ribo", assay = "RNA")

print("Plot continuous")

continuous_vars <- c("nFeature_RNA", "nCount_RNA", "pH", "PMI", "Age", "percent.mt", "percent.ribo")
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/Finalized_outputs/5_UMAPS_Harmony_Seeded_continuous.pdf",
 onefile = TRUE, width = 12, height = 16)
FeaturePlot(harmonized_object, reduction= "UMAPHarmonySeededBatchSampleChemistry",
                                               ncol = 2, features = continuous_vars)
dev.off()


#Cluster using optimal parameters from scclusteval

harmonized_object <- FindNeighbors(harmonized_object, reduction = "Harmony_Seeded_Batch_Sample_Chemistry", dims = 1:70, k.param = 30,
        nn.eps = 0, verbose = TRUE, force.recalc = TRUE)
harmonized_object <- FindClusters(object = harmonized_object, n.start = 100, resolution = 0.7, verbose = TRUE, algorithm = 2)


#Compare to scclusteval stable results

library(dplyr)
library(mclust)
library(scclusteval)
library(pheatmap)

fullsample_idents <- readRDS(file =
  "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/gather_full_sample.rds")
fullsample_idents %>%
  mutate(cluster_num = purrr::map_dbl(original_ident_full, ~n_distinct(.x))) -> fullsample_idents

#ARIs for all
print("Adjusted random indices for comparison")
fullsample_idents %>% mutate(ARI_seeded = map_dbl(original_ident_full, ~adjustedRandIndex(unlist(.x), Idents(harmonized_object)))) %>%
        select(pc, resolution, k_param, cluster_num, ARI_seeded)-> ARI_fullsample_ident
write.csv(ARI_fullsample_ident, 
file = "/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/Finalized_outputs/5_Seeded_Harmony_ARI_comparison.csv")

#PairWiseJaccardSets for identical parameters
fullsample_idents_selected <- fullsample_idents %>% filter(resolution == 0.7, k_param == 30, pc == 70) %>%
        select(original_ident_full)
JaccardSets_matrix <- PairWiseJaccardSets(unlist(fullsample_idents_selected$original_ident_full), Idents(harmonized_object))
pdf(file = 
"/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/Finalized_outputs/5_JaccardSets_Seeded_Unseeded_Optimal_Clustering.pdf", onefile = TRUE)
pheatmap(JaccardSets_matrix, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

saveRDS(harmonized_object, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_merging_normalization/5_reharmonized_seeded.Rds")

print("Print warnings summary.")
summary(warnings())

print("Print session info.")
sessionInfo()

