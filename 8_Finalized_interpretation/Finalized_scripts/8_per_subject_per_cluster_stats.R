library(Seurat)
library(tidyverse)

harmonized_object <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

print("Add back BRETIGEA module score results")
bret_scores <- read_csv(file = 
	"/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_outputs/2_BRETIGEA_BrainInABlender_scores.csv")
bret_scores <- as.data.frame(bret_scores)
rownames(bret_scores) <- colnames(harmonized_object)

harmonized_object <- AddMetaData(harmonized_object, bret_scores)

cluster_corr <- read_csv(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Cluster_correspondence.csv")
rename_list <- cluster_corr %>% select(Cluster, `RNA_snn_res.0.7`) %>% deframe
harmonized_object$RNA_snn_res.0.7 <- recode_factor(harmonized_object$RNA_snn_res.0.7, !!!rename_list)

base_path <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/"

inhib_idents <- levels(Idents(harmonized_object)) %>% grep(x = ., pattern = "^In", value = TRUE)
ol_idents <- c(levels(Idents(harmonized_object)) %>% grep(x = ., pattern = "^Ast", value = TRUE), 
	levels(Idents(harmonized_object)) %>% grep(x = ., pattern = "^O", value = TRUE))

print("Remake BRETIGEA, QC, and marker VlnPlots and also save the corresponding data")
pdf(file = paste0(base_path, "8_BRETIGEA_vlns_remade.pdf"), 
	height = 16, width = 14) 
VlnPlot(harmonized_object, features = c("ast", "end", "mic", "neu", "oli", "opc"), pt.siz = -1, 
	ncol =2, group.by = "RNA_snn_res.0.7") * theme(axis.text.x = element_text(angle = 90))
p1 <- FetchData(harmonized_object, vars = c("ast", "end", "mic", "neu", "oli", "opc", "RNA_snn_res.0.7")) 
write.csv(p1, file = paste0(base_path, "8_source_data_BRETIGEA_remade.csv"))

print("Inhibitory violins and source data")
VlnPlot(harmonized_object, features = c("VIP", "LAMP5", "GAD1", "GAD2", "LHX6", "ADARB2", "SST", "PVALB"),
	ncol = 2, pt.size = -1, idents = inhib_idents) * theme(axis.text.x = element_text(angle = 90)) 
p1 <- FetchData(harmonized_object, vars = c("VIP", "LAMP5", "GAD1", "GAD2", "LHX6", "ADARB2", "SST", "PVALB", "Cluster"),
        cell = WhichCells(harmonized_object, ident = inhib_idents))
write.csv(p1, file = paste0(base_path, "8_source_data_inhib_vlns_remade.csv"))

print("Oligo violins and source data")
VlnPlot(harmonized_object, features = c("PDGFRA", "ZFPM2", "ITPR2", "OLIG1", "MOG", "TCF7L2"), 
	ncol =2, pt.size = -1, idents = c(ol_idents, "End1", "Mic1")) * theme(axis.text.x = element_text(angle = 90)) 
VlnPlot(harmonized_object, features = c("PDGFRA", "ZFPM2", "ITPR2", "OLIG1", "MOG", "TCF7L2"),
        ncol =2, pt.size = -1, idents = ol_idents) * theme(axis.text.x = element_text(angle = 90))
p1 <- FetchData(harmonized_object, vars = c("PDGFRA", "ZFPM2", "ITPR2", "OLIG1", "MOG", "TCF7L2", "Cluster"),
	cells = WhichCells(harmonized_object, ident = ol_idents))
write.csv(p1, file = paste0(base_path, "8_source_data_ol_vlns_remade.csv"))

print("QC violins and source data")
VlnPlot(harmonized_object, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), log = TRUE,
        ncol =1, pt.size = -1, group.by = "RNA_snn_res.0.7") * theme(axis.text.x = element_text(angle = 90))
p1 <- FetchData(harmonized_object, vars = c("nCount_RNA", "nFeature_RNA", "percent.mt", "RNA_snn_res.0.7"))
write.csv(p1, file = paste0(base_path, "8_source_data_QC_vlns_remade.csv")) 
dev.off()

print("Generate tables per subject and cell type")
harmonized_object@meta.data %>% group_by(Sample) %>% 
		summarise(MedianUMI = median(nCount_RNA), MedianGene = median(nFeature_RNA), NumCells = n(), 
			Sex = unique(Sex), Condition = unique(Condition)) -> subject_stats
harmonized_object@meta.data %>% group_by(Cluster) %>% 
		 summarise(MedianUMI = median(nCount_RNA), MedianGene = median(nFeature_RNA), NumCells = n()) -> cluster_stats
harmonized_object@meta.data %>% group_by(Broad) %>%
                 summarise(MedianUMI = median(nCount_RNA), MedianGene = median(nFeature_RNA), NumCells = n()) -> broad_stats

print("Save UMAP coords as source data")
umap_coords <- Embeddings(harmonized_object, reduction="UMAPHarmonySeededBatchSampleChemistry")
umap_coords %>% as_tibble(rownames = "Cell") -> umap_coords
meta_to_add <- harmonized_object@meta.data[, c("Broad", "Cluster", "Sex", "Chemistry")]
meta_to_add %>% rownames_to_column(var = "Cell") -> meta_to_add
umap_coords <- inner_join(umap_coords, meta_to_add)

umap_coords_pca <- Embeddings(harmonized_object, reduction="UMAP_pca") %>% as_tibble(rownames = "Cell")
umap_coords <- inner_join(umap_coords, umap_coords_pca)

print("Generate more tables per group and sex")
subject_stats %>% group_by(Sex) %>% summarise(MedMedUMI = median(MedianUMI), MedMedGene = median(MedianGene))

base_path <- "/home/malosree/projects/def-gturecki/malosree/Finalized_interpretation_20220309/Finalized_outputs/"

write.csv(umap_coords, file = paste0(base_path, "8_source_data_umap_coords.csv"))

harmonized_object@meta.data %>% filter(Cluster %in%  c("InN1_PV", "InN9_PV", "Mic1") & Sex == "Female") %>% group_by(Cluster, Sample) %>%
                 summarise(NumCells = n()) -> female_validation_cluster_cell_num

write.csv(subject_stats, file = paste0(base_path, "8_subject_stats.csv")) 
write.csv(cluster_stats, file = paste0(base_path, "8_cluster_stats.csv"))
write.csv(broad_stats, file = paste0(base_path, "8_broad_stats.csv"))
write.csv(female_validation_cluster_cell_num, file = paste0(base_path, "8_female_validation_cluster_cell_num.csv"))

print("Make proportion plots")
Subject_props <- table(harmonized_object$Cluster, harmonized_object$Sample) %>% as.data.frame() %>% setNames(nm = c("Cluster", "Subject_props", "Freq"))
Batch_props <- table(harmonized_object$Cluster,  harmonized_object$Batch) %>% as.data.frame() %>% setNames(nm = c("Cluster", "Batch_props", "Freq"))
Chemistry_props <- table(harmonized_object$Cluster, harmonized_object$Chemistry) %>% as.data.frame() %>% setNames(nm = c("Cluster", "Chemistry_props", "Freq"))
Sex_props <- table(harmonized_object$Cluster, harmonized_object$Sex) %>% as.data.frame() %>% setNames(nm = c("Cluster", "Sex_props", "Freq"))
Condition_props <- table(harmonized_object$Cluster, harmonized_object$Condition) %>% as.data.frame() %>% setNames(nm = c("Cluster", "Condition_props", "Freq"))
Sex_condition_props <- table(harmonized_object$Cluster, paste(harmonized_object$Sex, harmonized_object$Condition)) %>% as.data.frame() %>% 
	setNames(nm = c("Cluster", "Sex_condition_props", "Freq"))
Male_condition_props <- table(harmonized_object$Cluster[harmonized_object$Sex == "Male"], 
		paste(harmonized_object$Sex, harmonized_object$Condition)[harmonized_object$Sex == "Male"]) %>% as.data.frame() %>%
		setNames(nm = c("Cluster", "Male_condition_props", "Freq")) %>% filter(Cluster != "ExN17")
Female_condition_props <- table(harmonized_object$Cluster[harmonized_object$Sex == "Female"],
                paste(harmonized_object$Sex, harmonized_object$Condition)[harmonized_object$Sex == "Female"]) %>% as.data.frame() %>%
                setNames(nm = c("Cluster", "Female_condition_props", "Freq"))
Bank_props <- table(harmonized_object$Cluster, harmonized_object$Bank) %>% as.data.frame() %>% setNames(nm = c("Cluster", "Bank_props", "Freq"))
Female_bank_props <- table(harmonized_object$Cluster[harmonized_object$Sex == "Female"],
                harmonized_object$Bank[harmonized_object$Sex == "Female"]) %>% as.data.frame() %>%
                setNames(nm = c("Cluster", "Female_bank_props", "Freq"))

pdf(file = paste0(base_path, "8_cluster_props.pdf"), onefile = TRUE, height = 5, width = 10)
	for(data_name in ls(pattern = "props$")) {
		data <- get(data_name)
		this_plot <- (ggplot(data, aes(fill=data[,data_name], y=Freq, x=Cluster)) + 
			geom_bar(position="fill", stat="identity") + labs(fill = data_name) +
			theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle = 90, size = 16)))
		print(this_plot)
		write.csv(this_plot$data, file = paste0(base_path, "8_source_data_", data_name, ".csv"))
	}
	Subject_props_norm <- table(harmonized_object$Cluster, harmonized_object$Sample)
	Subject_props_norm <- Subject_props_norm/rowSums(Subject_props_norm)  
	Subject_props_norm <- as.data.frame(Subject_props_norm) %>% setNames(n = c("Cluster", "Subject_props_norm", "Freq")) 
	this_plot <- (ggplot(Subject_props_norm, aes(y=Freq, x=Cluster, fill = Cluster)) +
                        geom_boxplot() +
                        theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle = 90, size = 16), legend.position = "none")) 
	print(this_plot)
	write.csv(this_plot$data, file = paste0(base_path, "8_source_data_subject_props.csv"))	
dev.off()
pdf(file = paste0(base_path, "8_batch_effect_UMAPS_redone.pdf"), onefile = TRUE, height = 7, width = 7)
	(DimPlot(harmonized_object, reduction = "UMAP_pca", group.by = "Sex", raster = FALSE) +
		ggplot2::theme_classic(base_size = 20)) %>% ggrastr::rasterize()
	(DimPlot(harmonized_object, reduction = "UMAP_pca", group.by = "Chemistry", raster = FALSE) +
		ggplot2::theme_classic(base_size = 20)) %>% ggrastr::rasterize()
	(DimPlot(harmonized_object, reduction = "UMAPHarmonySeededBatchSampleChemistry", group.by = "Sex", raster = FALSE) +
                ggplot2::theme_classic(base_size = 20)) %>% ggrastr::rasterize()
        (DimPlot(harmonized_object, reduction = "UMAPHarmonySeededBatchSampleChemistry", group.by = "Chemistry", raster = FALSE) +
                ggplot2::theme_classic(base_size = 20)) %>% ggrastr::rasterize()

dev.off()

print(summary(warnings()))
sessionInfo()

