library(Seurat)
library(ggplot2)
library(patchwork)

harmonized_object <- readRDS(file = "/home/malosree/scratch/Finalized_merging_normalization/3_harmonized_object.Rds")


#Plot how much UMAP corresponds to batch and other categorical variables.
print("Plot categoricals")

plots_list <- list()

reductions <- grep(names(harmonized_object@reductions), pattern = "^UMAP", value = TRUE)

categoricals  <- c("Batch", "Chemistry", "Sex", "Sample")

for(reduction in reductions){
	for(categorical in categoricals) {
		plots_list[[paste(reduction, categorical, sep = "_")]] <- DimPlot(harmonized_object, reduction= reduction, 
											group.by = categorical)
	}
}

pdf(file = "/home/malosree/scratch/Finalized_merging_normalization/Finalized_outputs/4_UMAPS_PCS_harmony_categorical.pdf",
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

plots_list <- list()

continuous_vars <- c("nFeature_RNA", "nCount_RNA", "pH", "PMI", "Age", "percent.mt", "percent.ribo")

for(reduction in reductions){
        plots_list[[reduction]] <- FeaturePlot(harmonized_object, reduction= reduction,
                                               ncol = 2, features = continuous_vars)
}

pdf(file = "/home/malosree/scratch/Finalized_merging_normalization/Finalized_outputs/4_UMAPS_PCS_harmony_continuous.pdf",
        onefile = TRUE, width = 12, height = 16)
plots_list
dev.off()

#Plot some basic QC parameters

print("Plot violins")
plots_list <- list()

features_to_plot <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")

log_scale <- c(TRUE, FALSE)

grouping_var <- c("Sample", "Chemistry", "Sex", "Batch")

for(var in grouping_var) {
	for(scale in log_scale) {
		plots_list[[paste(var, scale, sep = "_")]]<- VlnPlot(harmonized_object, features = features_to_plot, 
						log = scale, group.by = var,
						ncol = 2 , 
						pt.size=-1)*theme(axis.text.x = element_text(angle = 90))
	}
}

pdf(file = "/home/malosree/scratch/Finalized_merging_normalization/Finalized_outputs/4_nCount_nFeature_pt.mito_pt.ribo.pdf", 
	onefile = TRUE, height = 10, width = 12)
plots_list
dev.off()

#Plot QC parameters as scatter plots

print("Plot some scatter plots")

Idents(harmonized_object) <- harmonized_object$Batch
pdf(file = "/home/malosree/scratch/Finalized_merging_normalization/Finalized_outputs/4_Pt_mito_nCount_nFeature_scatter.pdf", 
	height = 12, width = 10)
plot1 <- FeatureScatter(harmonized_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(harmonized_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1/plot2
dev.off()

print("Print warnings summary")
summary(warnings())

print("Print session info.")
sessionInfo()

