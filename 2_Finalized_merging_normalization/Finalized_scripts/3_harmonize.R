library(Seurat)
library(harmony)

objects_merged <- readRDS(file = "/home/malosree/scratch/Finalized_merging_normalization/2_Finalized_merged_norm_obj.Rds")

print("Run PCA")
objects_merged <- RunPCA(objects_merged, npcs = 100)

#Using this approach based on this thread: https://github.com/immunogenomics/harmony/issues/41
#Was also in the (old) Signac tutorial here: https://satijalab.org/signac/0.2/articles/integration.html

#Batch correct PCA coordinates using Harmony, different combinations tested
#These technically return Seuratv3 objects. Checking version after running this step still says Seurat 4.0.1.
#Note that this will treat each run of S250 separately. 
print("Harmonizing")
variables_list <- list(c("Batch"), c("Chemistry"), c("Sample"), c("Batch", "Chemistry"), c("Batch", "Sample"), 
				c("Sample", "Chemistry"), c("Batch", "Sample", "Chemistry"))
pdf(file = "/home/malosree/scratch/Finalized_merging_normalization/Finalized_outputs/3_harmony_convergence_plots.pdf", onefile =TRUE)
for(items in variables_list) {
	name <- paste(c("Harmony", items), collapse = "_")
	objects_merged <- RunHarmony(objects_merged, reduction.save = name, group.by.vars = items, kmeans_init_iter_max = 5000, 
		plot_convergence = TRUE)
	print(name)
	print(summary(warnings()))
}
dev.off()

#just in case
harmonized_object <- UpdateSeuratObject(objects_merged)

#Elbow plots
print("Plotting elbows")

reductions <- names(harmonized_object@reductions)

plots_list <- list()

pdf(file = "/home/malosree/scratch/Finalized_merging_normalization/Finalized_outputs/3_Elbow_plots.pdf", onefile = TRUE)
for(reduction in reductions) {
	plots_list[[reduction]] <- ElbowPlot(harmonized_object, reduction = reduction, ndims = 100)
}
plots_list
dev.off()

#Run UMAPS using all calculated dimensions
print("Running UMAPS")
for(reduction in reductions){
	harmonized_object <- RunUMAP(harmonized_object, dims = 1:100, reduction= reduction, 
					reduction.name = paste("UMAP", reduction, sep = "_"),
					reduction.key = paste("UMAP", reduction, sep = "_"))
}

saveRDS(harmonized_object, file = "/home/malosree/scratch/Finalized_merging_normalization/3_harmonized_object.Rds")

#Check UMAP variation by categorial variable numerically
print("UMAP (all 100 dimensions) lms")
reductions <- grep(names(harmonized_object@reductions), pattern = "^UMAP", value = TRUE)
variables <- c("Batch", "Chemistry", "Sample")

var_results <- vector()

for(reduction in reductions) {
        cell.embeddings <- Embeddings(harmonized_object, reduction = reduction)
        print(reduction)
        results_1 <- vector()
        results_2 <- vector()
        for(variable in variables) {
                print(variable)
                results_1 <- c(results_1, summary(lm(cell.embeddings[,1]~harmonized_object@meta.data[, variable]))$r.squared)
                results_2 <- c(results_2, summary(lm(cell.embeddings[,2]~harmonized_object@meta.data[, variable]))$r.squared)
                if(length(results_2) == length(variables)) var_results <- rbind(var_results, results_1, results_2)
        }
}

rownames(var_results) <- paste(rep(reductions, each = 2), rep(1:2, times = length(reductions), sep = "_"))
colnames(var_results) <- variables

write.csv(var_results, file = "/home/malosree/scratch/Finalized_merging_normalization/Finalized_outputs/3_UMAP_categorical_vars_lm_all_100.csv")

print("Print warnings summary.")
summary(warnings())

print("Print session info.")
sessionInfo()

