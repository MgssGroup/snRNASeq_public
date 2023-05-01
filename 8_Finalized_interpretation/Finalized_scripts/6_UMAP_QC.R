library(Seurat)
library(tidyverse)

harmonized_object <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

print(table(harmonized_object$Cluster))

print("Get layer predictions from two sections in the Maynard data and add as metadata")
layer_preds <- read_csv(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/3.4_151673_Layer__preds.csv")
layer_preds %>% select(1:2) %>% column_to_rownames("...1") %>% as.data.frame() -> layer_preds

sum(rownames(layer_preds) == rownames(harmonized_object@meta.data))

harmonized_object <- AddMetaData(harmonized_object, metadata= layer_preds)

layer_preds <- read_csv(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/3.4_151507_Layer__preds.csv")
layer_preds %>% select(1:2) %>% column_to_rownames("...1") %>% as.data.frame() %>% setNames(n = c("predicted.id2"))-> layer_preds

sum(rownames(layer_preds) == rownames(harmonized_object@meta.data))

harmonized_object <- AddMetaData(harmonized_object, metadata= layer_preds)

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/6_additional_UMAPs.pdf", onefile = TRUE,
                height = 6, width = 8)
print("Remake broad and cluster UMAPs")
dp <- DimPlot(harmonized_object, reduction = "UMAPHarmonySeededBatchSampleChemistry", label = TRUE, raster = TRUE, group.by = "Broad", repel = TRUE, label.size = 6) + 
	ggplot2::theme_classic(base_size = 20) +  theme(legend.position = "none", axis.title = element_text(size = 0))
dp

#write.csv(dp$data, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/6_broad_umap_source_data.csv")

dp <- DimPlot(harmonized_object, reduction = "UMAPHarmonySeededBatchSampleChemistry", label = TRUE, raster = TRUE, 
	#group.by = "Cluster", 
	repel = TRUE, label.size = 5) +
        ggplot2::theme_classic(base_size = 20) + theme(legend.position = "none", axis.title = element_text(size = 0))
dp

#write.csv(dp$data, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/6_cluster_umap_source_data.csv")

harmonized_object@meta.data$Bank <- factor(harmonized_object@meta.data$Bank, levels = c("NSU", "Douglas"))
DimPlot(harmonized_object, reduction = "UMAPHarmonySeededBatchSampleChemistry", raster = FALSE,
        group.by = "Bank") +
        ggplot2::theme_classic(base_size = 20) + theme(axis.title = element_text(size = 0))
#harmonized_object <- subset(harmonized_object, Broad != "Mix" & Cluster != "ExN5" & Cluster != "ExN17")

#DimPlot(harmonized_object, reduction = "UMAPHarmonySeededBatchSampleChemistry", label = TRUE, raster = TRUE, group.by = "Broad") + ggplot2::theme(legend.position = "none")
print("Run UMAP variation")
harmonized_object <-  RunUMAP(harmonized_object, reduction = "Harmony_Seeded_Batch_Sample_Chemistry", dims = 1:50, reduction.name = "UMAP_50", n.neighbors = 50)

DimPlot(harmonized_object, reduction = "UMAP_50", label = TRUE, raster = TRUE, group.by = "Broad") + ggplot2::theme(legend.position = "none")
print("Run excitatory only UMAP")
harmonized_ex <- subset(harmonized_object, Broad == "ExN")

harmonized_ex <- RunUMAP(harmonized_ex, reduction = "Harmony_Seeded_Batch_Sample_Chemistry", dims = 1:50, reduction.name = "UMAP_Ex", n.neighbors = 20)
print("Run inhibitory only UMAP")
harmonized_in <- subset(harmonized_object, Broad == "InN")

harmonized_in <- RunUMAP(harmonized_in, reduction = "Harmony_Seeded_Batch_Sample_Chemistry", dims = 1:50, reduction.name = "UMAP_In", n.neighbors = 20)
print("Run glia only UMAP")
harmonized_glia <- subset(harmonized_object, Broad %in% c("Ast", "OPC", "Oli", "End", "Mic"))

harmonized_glia <- RunUMAP(harmonized_glia, reduction = "Harmony_Seeded_Batch_Sample_Chemistry", dims = 1:50, reduction.name = "UMAP_Glia", n.neighbors = 20)
print("Plot excitatory UMAPs labelled by layer and save the figure source data")
DimPlot(harmonized_ex, reduction = "UMAP_Ex", label = TRUE, raster = TRUE, repel = TRUE, label.size = 5) + theme_classic(base_size = 20) +  theme(legend.position = "none")
harmonized_ex <- SetIdent(harmonized_ex, value = "predicted.id")
DimPlot(harmonized_ex, reduction = "UMAP_Ex", label = FALSE, raster = TRUE, repel = TRUE, label.size = 5)  + theme_classic(base_size = 18) 
harmonized_ex <- SetIdent(harmonized_ex, value = "predicted.id2")
DimPlot(harmonized_ex, reduction = "UMAP_Ex", label = FALSE, raster = TRUE, repel = TRUE, label.size = 5)  + theme_classic(base_size = 18)

ex_UMAP <- Embeddings(harmonized_ex, reduction = "UMAP_Ex") %>% as_tibble(rownames = "Cell") 
meta_to_add <- harmonized_ex@meta.data[, c("Cluster", "predicted.id")] %>% as_tibble(rownames = "Cell")
ex_UMAP <- inner_join(ex_UMAP, meta_to_add)
write.csv(ex_UMAP, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/6_source_data_umap_ex.csv")
rm(ex_UMAP)
print("Plot the inhibitory UMAP labelled by layer and save the figure source data")
DimPlot(harmonized_in, reduction = "UMAP_In", label = TRUE, raster = TRUE, repel = TRUE, label.size = 5)  + theme_classic(base_size = 18) +  theme(legend.position = "none")

harmonized_in <- SetIdent(harmonized_in, value = "predicted.id2")
DimPlot(harmonized_in, reduction = "UMAP_In", label = FALSE, raster = TRUE, repel = TRUE, label.size = 5)  + theme_classic(base_size = 18)

in_UMAP <- Embeddings(harmonized_in, reduction = "UMAP_In") %>% as_tibble(rownames = "Cell")
meta_to_add <- harmonized_in@meta.data[, c("Cluster", "Broad")] %>% as_tibble(rownames = "Cell")
in_UMAP <- inner_join(in_UMAP, meta_to_add)
write.csv(in_UMAP, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/6_source_data_in_umap.csv")

DimPlot(harmonized_glia, reduction = "UMAP_Glia", label = TRUE, raster = TRUE)
FeaturePlot(harmonized_glia, features = c("PDGFRB", "CSF1R", "CLDN5", "VIM"), raster = TRUE, reduction = "UMAPHarmonySeededBatchSampleChemistry")

dev.off()
print("Write out the cluster tree")
ape::write.nexus(harmonized_object@tools$BuildClusterTree, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/6_dendrogram.nex")

print("Make DotPlot and write source data")
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/6_additional_dotplot.pdf", onefile = TRUE,
                height = 9, width = 11)
dp <- DotPlot(harmonized_object, features = c("SNAP25", "SLC17A7", "GAD1", "ALDH1L1", "PDGFRA", "PLP1", "CLDN5", "CX3CR1"), 
	ident = c("InN9_PV", "InN1_PV", "ExN17", "ExN12_L56", "ExN11_L56", "ExN13_L56", "ExN16_L56", "ExN20_L56", "ExN18", "ExN14", "ExN1_L24", "ExN3_L46", 
	"ExN6", "ExN2_L23", "ExN9_L23", "ExN8_L24", "ExN10_L46", "Ast1", "Ast2", 
	"InN10_ADARB2", "InN8_ADARB2", "InN3_VIP", "InN4_VIP", "ExN15_L56", "ExN19_L56", "InN5_SST", "InN2_SST", "InN6_LAMP5", 
	"InN7_Mix", "OPC3", "Oli3", "Oli1", "Oli2", "OPC1", "OPC2", "ExN7", "ExN4_L35", "ExN5", "Mix", "End1", "Mic1")) + theme_classic(base_size = 18)+
	theme(axis.text.x = element_text(angle = 90))
dp

write.csv(dp$data, file = "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/6_DotPlot_source_data.csv")

#the previous BuildClusterTree calls actually used the average expression of the the variable genes as dims was not specified, according to the documentation
harmonized_object <- BuildClusterTree(harmonized_object, dims = 1:70, reduction = "Harmony_Seeded_Batch_Sample_Chemistry", reorder = TRUE)

PlotClusterTree(harmonized_object)

dev.off()

print("Make subject proportions plots and save the source data")

cluster_cats <- c("Cluster", "Broad")

plots_list <- list()

harmonized_all <- harmonized_object

for (cluster_cat in cluster_cats) {

for(sex in list(c("Male", "Female"), c("Male"), c("Female"))) {

harmonized_object <- subset(harmonized_all, Sex %in% sex)

#Get subject proportions in each cluster
print("Get subject proportions in each cluster")
subject_proportions <- table(harmonized_object$Sample, t(harmonized_object@meta.data[cluster_cat]))

if("Male" %in% sex) {
	subject_proportions["M24",] <- subject_proportions["M24",]+subject_proportions["M24_2",]
	subject_proportions <- subject_proportions[setdiff(rownames(subject_proportions), "M24_2"),]
}

subject_proportions <- subject_proportions/rowSums(subject_proportions)

print(rowSums(subject_proportions))

subject_proportions %>% as_tibble(.name_repair = "minimal") %>% setNames(c("Subject", "Cluster", "Proportion")) %>%
                pivot_wider(names_from = Cluster, values_from = Proportion) %>%
                mutate(Condition = map_chr(Subject, ~unique(harmonized_object$Condition[harmonized_object$Sample == .x]))) -> subject_proportions

sex <- paste(sex, collapse = "_")

write.csv(subject_proportions, file = paste0(
	"/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/6_wilcoxon_source_data",
	sex, "_", cluster_cat,
	".csv"))

for(name in names(select_if(subject_proportions, is.numeric))) {
        print(name)
        myData <- data.frame(Proportion = subject_proportions[name][[1]], Condition = as.factor(subject_proportions$Condition),
                        Subject = subject_proportions$Subject)
        res_plot <- ggplot(data = myData, aes(x= Proportion, y = Condition, label = Subject, color = Condition))+geom_boxplot()+geom_text(size = 6)+coord_flip()+
                        ggtitle(name)+
                        theme_classic(base_size = 20)
        plots_list[[paste(name, sex, sep = "_")]] <- res_plot

	}

}
}

#Print boxplots out to file
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/6_boxplots_wilcoxon.pdf",
        onefile = TRUE)

        plots_list

dev.off()

print(table(harmonized_all$Cluster))


print(summary(warnings()))
sessionInfo()

