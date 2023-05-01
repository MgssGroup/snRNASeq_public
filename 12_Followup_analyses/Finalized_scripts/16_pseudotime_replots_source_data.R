library(tidyverse)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)

load(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/2_pseudotime_fitted.RData")
reducedDims(harmonized_object)$OL_UMAP %>% as_tibble(rownames = "Cell") -> OL_UMAP
as.data.frame(colData(harmonized_object)[, c("Cluster", "slingPseudotime_1")]) %>% 
	as_tibble(rownames = "Cell") %>% inner_join(OL_UMAP) -> OL_UMAP
write.csv(OL_UMAP, 
	file = "/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_outputs/16_source_data_PT_UMAP.csv")
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_outputs/16_PT_plots_regen.pdf", 
	onefile = TRUE, height = 8, width = 12)
for(gene in c("PLP1", "TCF7L2", "PDGFRA")) {
	for(object_name in c("male", "female")) {
		sce <- get(object_name)
		p1 <- ggrastr::rasterize(plotSmoothers(sce, counts(sce), gene)+ggtitle(gene))
		print(p1)
		write.csv(p1$data, 
			file = paste0("/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_outputs/", 
				"16_source_data_smoother_PT_", gene, "_", object_name, ".csv")) 
	}
}
dev.off()

print(summary(warnings()))
sessionInfo()

