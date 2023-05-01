library(Seurat)
library(tidyverse)

harmonized_object <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

new_cell_names <- paste(colnames(harmonized_object), harmonized_object$Broad, harmonized_object$Cluster, sep = ".")

harmonized_object <- RenameCells(harmonized_object, old.names = colnames(harmonized_object), new.names = new_cell_names)

library(Matrix)
library(tidyverse)
data <- GetAssayData(harmonized_object, assay= "RNA", slot = "counts")

writeMM(data, file = "/home/malosree/scratch/GEO_upload_female_snRNAseq/combined_counts_matrix.mtx")
write.csv(colnames(data), file = "/home/malosree/scratch/GEO_upload_female_snRNAseq/combined_counts_matrix_cells_columns.csv", row.names = FALSE)
write.csv(rownames(data), file = "/home/malosree/scratch/GEO_upload_female_snRNAseq/combined_counts_matrix_genes_rows.csv", row.names = FALSE)

ucsc_folder <- "/home/malosree/projects/rrg-gturecki/MGSS_MDD_dlPFC_snRNAseq_UCSCCellBrowser/"
metadata_ucsc <- harmonized_object@meta.data %>% as_tibble(rownames = "Cell") %>% select(Cell, Cluster, Broad, Sex, Batch, Chemistry, Condition)
umap_ucsc <- Embeddings(harmonized_object, reduction = "UMAPHarmonySeededBatchSampleChemistry") %>% as_tibble(rownames = "Cell")
exprs_ucsc <- GetAssayData(harmonized_object, assay = "RNA", slot = "counts")
exprs_ucsc_norm <- GetAssayData(harmonized_object, assay = "RNA", slot = "data")

#dim(exprs_ucsc)
#exprs_ucsc_filt <- exprs_ucsc[rowSums(exprs_ucsc) > 0,]
#dim(exprs_ucsc_filt)

write.table(metadata_ucsc, file = paste0(ucsc_folder, "metadata.tsv"), row.names = FALSE)
write.table(umap_ucsc, file = paste0(ucsc_folder, "umap.coords.tsv"), row.names = FALSE)
saveRDS(exprs_ucsc, file = paste0(ucsc_folder, "exprMatrix_sparse.Rds"))
saveRDS(exprs_ucsc_norm, file = paste0(ucsc_folder, "exprMatrix_sparse_norm.Rds"))
#write.table(exprs_ucsc_filt, file = paste0(ucsc_folder, "exprMatrix_sparse.tsv"))

print(summary(warnings()))
sessionInfo()
