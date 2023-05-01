#using the following vignettes as guides
#https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html#overview-1
#https://statomics.github.io/tradeSeq/articles/tradeSeq.html
#https://statomics.github.io/tradeSeq/articles/fitGAM.html

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(grDevices)
library(tidyverse)
library(ggrastr)
library(patchwork)
library(rstatix)
library(tradeSeq)

#Load, subset, and convert data
print("Load data")
harmonized_object <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

#print("Subsetting data, excluding OPC# as it is a small cluster, < 100 cells, and clusters away from remaining cells in UMAP")
OL_clusters <- c("Oli1", "Oli2", "Oli3", "OPC1", "OPC2", "OPC3")
harmonized_object <- subset(harmonized_object, Cluster %in% OL_clusters)

#Re-run UMAP
print("Re-running UMAP with Harmony only on OL cells")
harmonized_object <- RunUMAP(harmonized_object, reduction = "Harmony_Seeded_Batch_Sample_Chemistry", dims = 1:10, reduction.name = "OL_UMAP", 
	min.dist=0.1, spread = 5, n.neighbors = 100)

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Pseudotime/2_OL_UMAP.pdf", onefile = TRUE)
	DimPlot(harmonized_object, reduction = "UMAPHarmonySeededBatchSampleChemistry")	%>% rasterize()
	DimPlot(harmonized_object, reduction = "OL_UMAP") %>% rasterize()
	DimPlot(harmonized_object, reduction = "OL_UMAP", group.by = "Condition") %>% rasterize()
	DimPlot(harmonized_object, reduction = "OL_UMAP", group.by = "Chemistry") %>% rasterize()
	DimPlot(harmonized_object, reduction = "OL_UMAP", group.by = "Sex") %>% rasterize()
	DimPlot(harmonized_object, reduction = "OL_UMAP", group.by = "Batch") %>% rasterize()
	DimPlot(harmonized_object, reduction = "OL_UMAP", group.by = "Sample") %>% rasterize()
	DimPlot(harmonized_object, reduction = "OL_UMAP", cells = WhichCells(harmonized_object, idents = "OPC2")) %>% rasterize()
	DimPlot(harmonized_object, reduction = "OL_UMAP", cells = WhichCells(harmonized_object, idents = "Oli3")) %>% rasterize()
dev.off()

male <- subset(harmonized_object, subset = Sex == "Male")
female <- subset(harmonized_object, subset = Sex == "Female")

genes <- c("SOX10", "PDGFRA", "MYT1", "ZFPM2", "OLIG1", "OLIG2")

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Pseudotime/2_OL_UMAP_genes.pdf", onefile = TRUE)
	FeaturePlot(harmonized_object, reduction  = "OL_UMAP", features = genes) 
	FeaturePlot(male, reduction = "OL_UMAP", features = genes) 
	FeaturePlot(female, reduction = "OL_UMAP", features = genes) 
dev.off()

genes2 <- c( "MBP", "PLP1", "MAG", "MOG", "TCF7L2", "GAPDH")

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Pseudotime/2_OL_UMAP_genes_2.pdf", onefile = TRUE)
        FeaturePlot(harmonized_object, reduction  = "OL_UMAP", features = genes2)
        FeaturePlot(male, reduction = "OL_UMAP", features = genes2)
        FeaturePlot(female, reduction = "OL_UMAP", features = genes2)
dev.off()

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Pseudotime/2_OL_UMAP_technical.pdf", onefile = TRUE)
        FeaturePlot(harmonized_object, reduction  = "OL_UMAP", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
        FeaturePlot(male, reduction = "OL_UMAP", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
        FeaturePlot(female, reduction = "OL_UMAP", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
dev.off()

#save.image(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/2_pseudotime.RData")

print("Convert to SingleCellExperiment")
harmonized_object <- as.SingleCellExperiment(harmonized_object)

harmonized_object <- harmonized_object[rowSums(counts(harmonized_object)>0)>0, ]

#perform slingshot trajectory analysis
print("Perform slingshot trajectory analysis")
harmonized_object <- slingshot(harmonized_object, reducedDim = "OL_UMAP", extend = "n", 
				clusterLabels = colData(harmonized_object)$Cluster, 
				start.clus = "OPC2", end.clus = "Oli3", 
				stretch = 0.1, 
				thresh = 0.3)

#Make dimensionality reduction plots
print("Make dimensionality reduction plots with pseudotime")

colors <- colorRampPalette(RColorBrewer::brewer.pal(9,'Spectral')[-6])(100)
plotcol <- colors[cut(harmonized_object$slingPseudotime_1, breaks = 100)]

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Pseudotime/2_slingshot_UMAP.pdf", onefile = TRUE)
	plot(reducedDims(harmonized_object)$OL_UMAP, col = plotcol, pch=16, asp =1)
	lines(SlingshotDataSet(harmonized_object), lwd=2, col='black')	
	plot(reducedDims(harmonized_object)$OL_UMAP, col = RColorBrewer::brewer.pal(6,'Set1')[as.numeric(harmonized_object$Cluster)-29], pch=16, asp = 1)
	lines(SlingshotDataSet(harmonized_object), lwd=2, type = 'lineages', col = 'black', show.constraints = TRUE)
dev.off()

print("Subset male and female data for differential expression analysis")
male <- harmonized_object[,harmonized_object$Sex == "Male"]
female <- harmonized_object[,harmonized_object$Sex == "Female"]

save.image(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/2_pseudotime.RData")

#compare if there are differences in pseudotime between cases and controls
print("Use KS tests to check for differences in pseudotime between cases and controls")

plots_list <- list()
res_ks <- list()

for(object_name in c("male", "female", "harmonized_object")) {
	sce <- get(object_name)
	condition <- colData(sce)$Condition
	data <- data.frame(Lineage1 = slingPseudotime(sce),condition = condition)
	plots_list[[object_name]] <- ggplot(data = data, aes(x=Lineage1, colour = condition)) + geom_density()
	res_ks[[object_name]] <- ks.test(slingPseudotime(sce)[colData(sce)$Condition == "Case", 1],
		slingPseudotime(sce)[colData(sce)$Condition == "Control", 1])
}

res_ks

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Pseudotime/2_density_plots_pseudotime.pdf", onefile = TRUE)
plots_list
dev.off()

object_names <- c("male", "female", "harmonized_object")

#Detailed pseudotime difference comparison between cases and controls

cluster_cats <- c("Cluster", "Broad")

#Function for running Wilxocon tests before bootstrapping
wilcoxon_PT <- function(median_PT, conditions, cluster){
        myData <- data.frame(median_PT = median_PT, Condition = conditions)
        res_test <- wilcox_test(formula = median_PT ~ Condition, data = myData, alternative = "two.sided", paired = FALSE, detailed = TRUE)
        res_test <- cbind(res_test, Cluster = cluster$Cluster)
}

plot_box <- function(data, id){
        res_plot <-  ggplot(data, aes(x= median_PT, y = Condition, label = Subject, color = Condition))+
        geom_boxplot()+geom_text()+
        coord_flip()+
        ggtitle(paste(object, id$Cluster))+
        theme_classic()
        res_plot
}

PT_wilcoxon_results_compiled <- list()

gene_count_plots_all <- list()

genes <- c("SOX10", "PDGFRA", "MYT1", "OLIG1", "OLIG2", "TCF7L2", "MBP", "PLP1", "GAPDH")

for(object in object_names){
        sce <- get(object)
        colData(sce)$Sample <- str_replace(colData(sce)$Sample, pattern = fixed("M24_2"), replacement = "M24")
        for (cluster_cat in cluster_cats) {
                print("Summarise the pseudotime by subject and cluster")
                PT_data <- tibble(Pseudotime = slingPseudotime(colData(sce)$slingshot),
                                Condition = as.factor(colData(sce)$Condition),
                                Subject = colData(sce)$Sample,
                                Cluster = colData(sce)[cluster_cat][[1]])
                PT_data %>% group_by(Cluster, Subject) %>%
                        summarise(median_PT = median(Pseudotime), Condition = unique(Condition)) -> PT_data

                print("Run Wilcoxon tests on pseudotime")
                PT_data %>% group_map(~wilcoxon_PT(.x$median_PT, .x$Condition, .y)) %>% bind_rows() -> results_wilcoxon
                PT_wilcoxon_results_compiled[[paste(object, cluster_cat, sep = "_")]] <- results_wilcoxon

                print("Make boxplots of pseudotime")
                PT_data %>% group_map(~plot_box(.x, .y)) -> plots_list[[paste(object, cluster_cat, sep = "_")]]

        }

        #Make some plots of selected genes along the pseudotime UMAP")
        print("Make gene count plots on pseudotime UMAP")
        gene_count_plots <- list()
        for(gene in genes) {
                gene_count_plots[[gene]] <- rasterize(plotGeneCount(sce, gene = gene))
        }
        gene_count_plots_all[[object]] <- Reduce("+", gene_count_plots) + plot_layout(ncol = 3)
}

#Write out Wilcoxon results
print("Write out Wilcoxon test results")
write.csv(bind_rows(PT_wilcoxon_results_compiled, .id = "Comparison"),
        file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Pseudotime/2_PT_wilcoxon.csv")

#Print boxplots out to file
print("Save plots")
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Pseudotime/2_boxplots_PT.pdf",
        onefile = TRUE)
        unlist(plots_list, recursive = FALSE)
dev.off()

#Save plots of selected genes along the pseudotime UMAP
print("Plot some genes on pseudotime UMAP")
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Pseudotime/2_pseudotime_marker_gene_expression.pdf",
        height = 12, width = 24)
        lapply(gene_count_plots_all, function(x){
        x
})
dev.off()

print(summary(warnings()))
sessionInfo()
