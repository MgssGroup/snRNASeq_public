library(Seurat)
library(harmony)
library(cluster)
library(muscat)
library(SingleCellExperiment)
library(tidyverse)
library(limma)
library(ggpmisc)
library(UpSetR)


print("Set base path, load and subset data")
base_path <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/"

harmonized_object <- readRDS("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

microglial_object <- subset(harmonized_object, Cluster == "Mic1")

print("Make some QC plots for clusters of interest which have DEGs in the male or female analysis")

clusters <- c("Mic1", "ExN10_L46", "Ast1", "Ast2", "OPC1", "OPC2", "InN1_PV", "InN9_PV", "InN8_Mix", "InN3_VIP", "ExN9_L23")
cluster_broad <- c("Ast", "OPC")
variables <- c("Chemistry", "Condition")

pdf(file= paste0(base_path, "Finalized_outputs/Microglia_refiltered/4.2_QC_plots_conitnued.pdf"))
	harmonized_object <- SetIdent(harmonized_object, value = harmonized_object$Cluster)
	DimPlot(harmonized_object, reduction = "UMAPHarmonySeededBatchSampleChemistry", label = TRUE, raster = TRUE) + theme(legend.position = "none")
	for(cluster in clusters) {
		for(var in variables) {
			(DimPlot(harmonized_object, cells = WhichCells(harmonized_object, idents = cluster), group.by = var,  
				reduction = "UMAPHarmonySeededBatchSampleChemistry", raster = TRUE) + ggtitle(cluster)) %>% print()
		}
	}
	harmonized_object <- SetIdent(harmonized_object, value = harmonized_object$Broad)
	DimPlot(harmonized_object, reduction = "UMAPHarmonySeededBatchSampleChemistry", label = TRUE, raster = TRUE)
	for(cluster in cluster_broad) {
                for(var in variables) {
                        (DimPlot(harmonized_object, cells = WhichCells(harmonized_object, idents = cluster), group.by = var,  
				reduction = "UMAPHarmonySeededBatchSampleChemistry", raster = TRUE) + ggtitle(cluster))%>% print()
                }
        }
dev.off()

rm(harmonized_object)

print("Find variables features, scale data, run PCA, harmonize.")
microglial_object <- FindVariableFeatures(microglial_object)
microglial_object <- ScaleData(microglial_object, vars.to.regress = c("nCount_RNA", "percent.mt"), model.use = "poisson", use.umi = TRUE)
microglial_object <- RunPCA(microglial_object, npcs = 50)

pdf(file = paste0(base_path, "Finalized_outputs/Microglia_refiltered/4.2_microglia_harmony_convergece.pdf"))
	set.seed(13)
	microglial_object <- RunHarmony(microglial_object, reduction.save = "MicHarmony", group.by.vars = c("Batch", "Sample", "Chemistry"), kmeans_init_iter_max = 5000,
                plot_convergence = TRUE)
dev.off()

microglial_object <- UpdateSeuratObject(microglial_object)

#Elbow plot
print("Plotting elbow")
pdf(file = paste0(base_path, "Finalized_outputs/Microglia_refiltered/4.2_Mic_Seeded_Harmony_Elbow_plot.pdf"), onefile = TRUE)
ElbowPlot(microglial_object, reduction = "MicHarmony", ndims = 50)
dev.off()

#Run UMAP
print("Running UMAP")
microglial_object <- RunUMAP(microglial_object, dims = 1:20, reduction= "MicHarmony",
                                        reduction.name = "UMAPMicHarmony",
                                        reduction.key = "UMAPMicHarmony")

zero_genes <- rownames(microglial_object)[rowSums(microglial_object) == 0]
length(zero_genes)
length(intersect(zero_genes, microglial_object@assays$RNA@var.features))
length(microglial_object@assays$RNA@var.features)

print("Cluster using a range of resolutions and dimensions")

clusterings_list <- list()

for(res in c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4)) {
        for(n_dims in c(15, 20, 25)) {
                microglial_object <- FindNeighbors(microglial_object, reduction = "MicHarmony", dims = 1:n_dims, verbose = FALSE)
                microglial_object <- FindClusters(microglial_object, resolution = res, verbose = FALSE)
		print(paste("params", res, n_dims, sep = "_"))
		table(microglial_object@meta.data[[paste0("RNA_snn_res.", res)]]) %>% print()
                clusterings_list[[paste("params", res, n_dims, sep = "_")]] <- as.integer(microglial_object@meta.data[[paste0("RNA_snn_res.", res)]])
        }
}

print("Get silhouette assessment for each clustering")
dist_mic <- dist(Embeddings(microglial_object, reduction = "MicHarmony"))

sils <- vector()

pdf(file = paste0(base_path, "Finalized_outputs/Microglia_refiltered/4.2_silhouettes_redone.pdf"), onefile = TRUE)
for(clustering in names(clusterings_list)) {
        this_sil <- silhouette(clusterings_list[[clustering]], dist_mic, FUN = median)
        #print(summary(this_sil)$avg.width)
        plot(this_sil, nmax.lab = 40, max.strlen = 5,
                main = clustering, sub = NULL, xlab = expression("Silhouette width "* s[i]),
                col = "gray",  do.col.sort = length(col) > 1, border = 0,
                cex.names = par("cex.axis"), do.n.k = TRUE, do.clus.stat = TRUE)
        sils <- rbind(sils, data.frame(Params = clustering, 
				Clus_Avg = summary(this_sil)$clus.avg.widths %>% mean(), 
				Avg = summary(this_sil)$avg.width))
}
ggplot(sils, aes(x= Avg, y = Clus_Avg, label = Params))+geom_label()+theme_classic()
dev.off()

max(sils$Clus_Avg)
sils$Params[sils$Clus_Avg== max(sils$Clus_Avg)]

max(sils$Avg)
sils$Params[sils$Avg== max(sils$Avg)]

print("Cluster with selected parameters")
microglial_object <- FindNeighbors(microglial_object, reduction = "MicHarmony", dims = 1:25)
microglial_object <- FindClusters(microglial_object, resolution = 0.01)

OL_lineage_markers <- c("MBP", "PLP1", "MAG", "MOG", "MOBP", "SOX10", "OLIG1", "OLIG2", "PDGFRA", "PCDH15", "TCF7L2", "MYT1", "ZFPM2")
Microglia_markers <- c("MRC1", "CX3CR1", "TMEM119", "FKBP5", "CXCR4", "P2RY12", "AIF1", "ITPR2")
Neuronal_markers <- c("SNAP25", "SLC17A7", "GAD1")	
plotting_features3 <- c("pANN", "nCount_RNA", "nFeature_RNA")

print("Make some plots of marker genes for the subclusters")
pdf(file = paste0(base_path, "Finalized_outputs/Microglia_refiltered/4.2_Redone_OL_Micro_markers_reclustered.pdf"), onefile = TRUE)
                DimPlot(microglial_object, reduction = "UMAPHarmonySeededBatchSampleChemistry", raster = TRUE)
                VlnPlot(microglial_object, features = Microglia_markers)
                VlnPlot(microglial_object, features = OL_lineage_markers)
		VlnPlot(microglial_object, features = Neuronal_markers)
		VlnPlot(microglial_object, features = plotting_features3)
		DimPlot(microglial_object, group.by = "Condition", reduction = "UMAPMicHarmony", raster = TRUE, split.by = "Sex")
		DimPlot(microglial_object, group.by = "Chemistry", reduction = "UMAPMicHarmony", raster = TRUE, split.by = "Sex")
		
dev.off()

table(Idents(microglial_object), microglial_object$Sex, microglial_object$Condition)

print("Get subcluster markers")

all_markers_microglia_presto <- presto::wilcoxauc(microglial_object, "RNA_snn_res.0.01")
all_markers_microglia_presto <- filter(as_tibble(all_markers_microglia_presto), padj < 0.05 & logFC > log(1.5) & pct_in-pct_out > 10)

write_csv(all_markers_microglia_presto, file = paste0(base_path, "Finalized_outputs/Microglia_refiltered/4.2_all_markers_microglia_Presto_WilcoxAUC.csv"))

write.csv(microglial_object@meta.data[,"RNA_snn_res.0.01"], file = paste0(base_path, "Finalized_outputs/Microglia_refiltered/4.2_microglial_thorough_subclusters.csv"))

print("Rerun differential expression using only cluster zero")

microglia_female_subset <- subset(microglial_object, idents = c("0"), subset = Sex == "Female")

table(Idents(microglia_female_subset), microglia_female_subset$Sex, microglia_female_subset$Condition)

directories <- c("Broad" = "03_Female_dataset_analysis/01_pseudobulk/01_cell_Broad/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_type/",
                 "Cluster" = "03_Female_dataset_analysis/01_pseudobulk/02_cell_subtype/10_Female_Pseudobulk_edgeR__Batch_Age_PMI_pH_cell_subtype/")

num_cells <- c("Broad" = 10, "Cluster" = 5)

for(category in c("Cluster", "Broad")) {
        print("Convert to SCE and prep SCE for muscat")
        sce <- as.SingleCellExperiment(microglia_female_subset)

        sce <- sce[rowSums(counts(sce) > 0) > 0, ]

        sce <- prepSCE(
                sce,
                kid = category,
                sid = "Sample",
                gid = "Condition",
                drop = FALSE)

        print("Prep the design matrix and contrasts. Save design matrix to file.")
        temp <- colData(sce)
        rownames(temp) <- NULL
        anno <-
                unique.data.frame(temp[c("sample_id",
                            "group_id",
                            "Sex",
                            "Batch",
                            "Age",
                            "pH",
                            "PMI",
                            "Chemistry")])
        anno$pH[anno$sample_id == "F35"] <- 6.47

        design <- model.matrix(~ 0+group_id+Batch+Age+PMI+pH,  data = anno)

        dimnames(design) <-
                list(anno$sample_id, levels = make.names(colnames(design)))

        write.csv(design, file = paste0(base_path, "/Finalized_outputs/Microglia_refiltered/",
                "4.2_subsetted_microglia_design_matrix_input.csv"))
        contrast <-
                makeContrasts(
                Case_vs_Control = (group_idCase - group_idControl),
                        levels = make.names(colnames(design))
                 )
	write.csv(contrast, file = paste0(base_path, "/Finalized_outputs/Microglia_refiltered/",
                "4.2_subsetted_microglia_contrast.csv"))
        print("Pseudobulk")
        pb <- aggregateData(
                sce,
                assay = "counts",
                fun = "sum",
                by = c("cluster_id", "sample_id"))

        print("Run differential expression analysis")
        res <- pbDS(
                pb,
                design = design,
                min_cells = num_cells[category],
                contrast = contrast,
                filter = "none",
                method = "edgeR")

        print("Format results. Write outputs.")
        results <- resDS(sce,
                res,
                cpm = TRUE,
                frq = FALSE,
                bind = "row",
                digits = 10)

        metadata <- read_csv(file = paste0(base_path, "Male_female_metadata_combined.csv")) %>%
                mutate(Sample = paste0(Sample, ".cpm"))

        case_subjects <- metadata %>% filter(Condition == "Case") %>% select(Sample) %>% unlist()
        print(case_subjects)
        control_subjects <- metadata %>% filter(Condition == "Control") %>% select(Sample) %>% unlist()
        print(control_subjects)

        results %>% mutate(NumNonZero = rowSums(across(ends_with(".cpm"), ~(.x>0)), na.rm = TRUE),
                        NumCaseExcluded = rowSums(across(any_of(case_subjects), is.na)),
                        NumControlExcluded = rowSums(across(any_of(control_subjects), is.na))) %>%
                        mutate(Greater3 = NumNonZero > 2) %>% mutate(score = logFC * (-log10(p_val))) -> results

        results %>% filter(p_adj.loc < 0.05, abs(logFC) > log2(1.1), Greater3 == TRUE) -> filtered_results

        write.csv(results, file = paste0(base_path, "Finalized_outputs/Microglia_refiltered/",
                                "4.2_Table_raw_masterfile_cpm_Freq_DS_stats_greater3_", category,".csv"))
        write.csv(filtered_results, file = paste0(base_path, "Finalized_outputs/Microglia_refiltered/",
                                "4.2_Table_raw_masterfile_cpm_Freq_DS_stats_FDR_05Pct_", category,".csv"))

        print(lapply(res$fit, function(x){sum(x$failed == TRUE)}))

        print("Compare to previous results")

        directory <- directories[category]

        results_prev <- read_csv(file = paste0(base_path, directory,
                                "01_Table_raw_masterfile_cpm_Freq_DS_stats_greater3.csv")) %>% filter(grepl("Mic", cluster_id)) %>%
                                mutate(score = logFC * (-log10(p_val)))

        results_prev %>% filter(p_adj.loc < 0.05, abs(logFC) > log2(1.1), Greater3 == TRUE) -> filtered_results_prev

        results %>% select(gene, logFC, score, p_val) -> results_for_corr
        results_prev %>% select(gene, logFC, score, p_val) -> results_prev
        this_df <- inner_join(results_for_corr, results_prev, by = "gene", suffix = c(".filtered", ".all"))
	this_df_selected <- this_df %>% arrange(p_val.all) %>% slice_head(n=3000)

	corr_plotting_func <- function(df) {
		this_df <- df
		this_res <- lm(logFC.filtered~logFC.all, data = drop_na(this_df)) %>% summary()
                this_annot <- list(R_squared = this_res$r.squared, p_val = this_res$coefficients[2,4], df = this_res$df[2]) %>% as_tibble()
                p1 <- (ggplot(data = drop_na(this_df), aes(x = logFC.all, y = logFC.filtered)) +
                        geom_smooth(method = "lm", se=FALSE, color="black") +
                        geom_point() +
                        annotate("table", label = this_annot, x= 0.2, y = 3.5) +
                        theme_classic(base_size = 22) + theme(axis.text = element_text(size = 22))) %>% ggrastr::rasterize() 
		print(p1)
		write.csv(p1$data, file = paste0(base_path,
			"Finalized_outputs/Microglia_refiltered/4.2_source_data_female_microglia_filter_no_filter_corrs", category,".csv"))
                this_res <- lm(score.filtered~score.all, data = drop_na(this_df)) %>% summary()
                this_annot <- list(R_squared = this_res$r.squared, p_val = this_res$coefficients[2,4], df = this_res$df[2]) %>% as_tibble()
                (ggplot(data = drop_na(this_df), aes(x = score.all, y = score.filtered)) +
                        geom_smooth(method = "lm", se=FALSE, color="black") +
                        geom_point() +
                        annotate("table", label = this_annot, x= 1, y = 30) +
                        theme_classic(base_size = 22) + theme(axis.text = element_text(size = 22))) %>% ggrastr::rasterize() %>% print()
	}
	
        print("Plot the correlations of logFCs and F stats with the previous results")
        pdf(file = paste0(base_path,
                "Finalized_outputs/Microglia_refiltered/4.2_female_microglia_filter_no_filter_corrs", category,".pdf"),
                onefile = TRUE, height = 10, width = 10)
		corr_plotting_func(this_df_selected)
		corr_plotting_func(this_df)
        dev.off()

        print("Print and plot overlaps with previous results at different thresholds")
        pdf(file = paste0(base_path,
                "Finalized_outputs/Microglia_refiltered/4.2_female_microglia_filter_no_filter_overlaps", category,".pdf"), onefile = TRUE, height = 10, width = 10)

        print(category)

        print("Cutoff 0.05")
        results %>% filter(p_adj.loc < 0.05, abs(logFC) > log2(1.1), Greater3 == TRUE) -> filtered_results_0.05
        intersect(filtered_results_prev$gene, filtered_results_0.05$gene) %>% print()
        length(filtered_results_0.05$gene) %>% print()

        print("Cutoff 0.1")
        results %>% filter(p_adj.loc < 0.1, abs(logFC) > log2(1.1), Greater3 == TRUE) -> filtered_results_0.1
        intersect(filtered_results_prev$gene, filtered_results_0.1$gene) %>% print()
        length(filtered_results_0.1$gene) %>% print()

        print("Cutoff 0.2")
        results %>% filter(p_adj.loc < 0.2, abs(logFC) > log2(1.1), Greater3 == TRUE) -> filtered_results_0.2
        intersect(filtered_results_prev$gene, filtered_results_0.2$gene) %>% print()
        length(filtered_results_0.2$gene) %>% print()

        upset(fromList(list(filtered_results_0.05 = filtered_results_0.05$gene,
                filtered_results_0.1 = filtered_results_0.1$gene,
                filtered_results_0.2 = filtered_results_0.2$gene,
                filtered_results_ori = filtered_results_prev$gene))) %>% print()

        dev.off()

}

print(summary(warnings()))
sessionInfo()
