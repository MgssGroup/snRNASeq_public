#load required libraries
library(Seurat)
library(dplyr)
library(readr)
library(stringr)

#list of male subjects, read in matrices
male_subjects <- setdiff(c(paste0("M", 1:34), "M24_2"), "M25")
for(subject in male_subjects){
 temp <- Read10X(data.dir = paste0("/home/malosree/scratch/male_counts/", subject, "_count/outs/filtered_feature_bc_matrix/"))
 assign(paste0("matrix_", subject), temp)
}

#list of female subjects, read in matrices
female_subjects <- paste0("F", 1:38)
for(subject in female_subjects){
 temp <- Read10X(data.dir = paste0("/home/malosree/scratch/female_counts/", subject, "_count/outs/filtered_feature_bc_matrix/"))
 assign(paste0("matrix_", subject), temp)
}

#read in metadata, check that data types make sense
meta_data <- read_csv(file = "/home/malosree/scratch/Male_female_metadata_combined.csv")
print("meta_data specs")
spec(meta_data)

#Create Seurat object, add metadata, and rename cells 
#technically some of these metadata are only for the aggregated metrics for the subject produced by CellRanger, 
#and not the set of cells we retain after filtering. 

for(subject in c(male_subjects, female_subjects)) {
 temp <- get(paste0("matrix_",subject))
 temp_meta <- filter(meta_data, Sample == subject)
 temp_meta <- slice(temp_meta, rep(1:n(), each = dim(temp)[2]))
 temp_meta <- as.data.frame(temp_meta)
 rownames(temp_meta) <- colnames(temp)
 temp_seurat <- CreateSeuratObject(counts = temp, meta.data = temp_meta)
 temp_seurat <- RenameCells(temp_seurat, new.names = paste(subject, colnames(temp), sep = fixed("."))) 
 assign(paste0("Seurat_", subject), temp_seurat)
 rm(list = c(paste0("matrix_",subject)))
}

#select female subjects by chemistry
meta_data %>% filter(Chemistry == "v2" & Sex == "Female") %>% select(Sample) -> female_subs_v2
meta_data %>% filter(Chemistry == "v3" & Sex == "Female") %>% select(Sample) -> female_subs_v3

male_stats_pre_filter <- data.frame()
female_v2_stats_pre_filter <- data.frame()
female_v3_stats_pre_filter <- data.frame()

#Calculate mitochondrial reads percentage, then remove mitochondrial genes, find the numbers of UMIs and genes per cell across all male subjects
print("Preprocessing male subjects")
for(subject in male_subjects) {
 temp_seurat <- get(paste0("Seurat_", subject))
 print(c(subject,dim(temp_seurat)))
 temp_seurat[['percent.mt']] <- PercentageFeatureSet(temp_seurat, pattern = "^MT-")
 genes_to_keep <- setdiff(rownames(temp_seurat), grep(rownames(temp_seurat), pattern = "^MT-", value = TRUE))
 temp_seurat <- subset(temp_seurat, features = genes_to_keep)
 male_stats_pre_filter <- rbind(male_stats_pre_filter, temp_seurat@meta.data[,c("Sample", "nFeature_RNA", "nCount_RNA")])
 assign(paste0("Seurat_",subject), temp_seurat)
}

#Calculate mitochondrial reads percentage, then remove mitochondrial genes, find the numbers of UMIs and genes per cell across v2 female subjects
print("Preprocessing female v2 subjects")
for(subject in female_subs_v2$Sample) {
 temp_seurat <- get(paste0("Seurat_", subject))
 print(c(subject,dim(temp_seurat)))
 temp_seurat[['percent.mt']] <- PercentageFeatureSet(temp_seurat, pattern = "^MT-")
 genes_to_keep <- setdiff(rownames(temp_seurat), grep(rownames(temp_seurat), pattern = "^MT-", value = TRUE))
 temp_seurat <- subset(temp_seurat, features = genes_to_keep)
 female_v2_stats_pre_filter <- rbind(female_v2_stats_pre_filter, temp_seurat@meta.data[,c("Sample", "nFeature_RNA", "nCount_RNA")])
 assign(paste0("Seurat_",subject), temp_seurat)
}

#Calculate mitochondrial reads percentage, then remove mitochondrial genes, find the numbers of UMIs and genes per cell across v3 female subjects
print("Preprocessing female v3 subjects")
for(subject in female_subs_v3$Sample) {
 temp_seurat <- get(paste0("Seurat_", subject))
 print(c(subject,dim(temp_seurat)))
 temp_seurat[['percent.mt']] <- PercentageFeatureSet(temp_seurat, pattern = "^MT-")
 genes_to_keep <- setdiff(rownames(temp_seurat), grep(rownames(temp_seurat), pattern = "^MT-", value = TRUE))
 temp_seurat <- subset(temp_seurat, features = genes_to_keep)
 female_v3_stats_pre_filter <- rbind(female_v3_stats_pre_filter, temp_seurat@meta.data[,c("Sample", "nFeature_RNA", "nCount_RNA")])
 assign(paste0("Seurat_",subject), temp_seurat)
}

#find quantiles for genes and UMIs dected per cell for each dataset and chemstry combination

print("male_stats_pre_filter nFeature_RNA")
quantile(male_stats_pre_filter$nFeature_RNA, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.98, 0.99, 0.995, 0.9975, 0.998, 0.999, 1))
print("male_stats_pre_filter nCount_RNA")
quantile(male_stats_pre_filter$nCount_RNA, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.98, 0.99, 0.995, 0.9975, 0.998, 0.999, 1))

print("female_v2_stats_pre_filter nFeature_RNA")
quantile(female_v2_stats_pre_filter$nFeature_RNA, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.98, 0.99, 0.995, 0.9975, 0.998, 0.999, 1))
print("female_v2_stats_pre_filter nCount_RNA")
quantile(female_v2_stats_pre_filter$nCount_RNA, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.98, 0.99, 0.995, 0.9975, 0.998, 0.999, 1))

print("female_v3_stats_pre_filter nFeature_RNA")
quantile(female_v3_stats_pre_filter$nFeature_RNA, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.98, 0.99, 0.995, 0.9975, 0.998, 0.999, 1))
print("female_v3_stats_pre_filter nCount_RNA")
quantile(female_v3_stats_pre_filter$nCount_RNA, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.98, 0.99, 0.995, 0.9975, 0.998, 0.999, 1))

male_stats_post_filter <- data.frame()
female_v2_stats_post_filter <- data.frame()
female_v3_stats_post_filter <- data.frame()

#Filter male datasets by gene, UMI count, and percent.mt, recalculate gene and UMI stats
print("Filtering male subjects")
for(subject in male_subjects) {
 temp_seurat <- get(paste0("Seurat_", subject))
 temp_seurat <- subset(temp_seurat, nCount_RNA < 35000 & nFeature_RNA > 350 & percent.mt < 10)
 print(c(subject,dim(temp_seurat)))
 male_stats_post_filter <- rbind(male_stats_post_filter, temp_seurat@meta.data[,c("Sample", "nFeature_RNA", "nCount_RNA")])
 assign(paste0("Seurat_",subject), temp_seurat)
}

#Filter v2 female datasets by gene, UMI counts, and percent.mt recalculate gene and UMI stats
print("Filtering female v2 subjects")
for(subject in female_subs_v2$Sample) {
 temp_seurat <- get(paste0("Seurat_", subject))
 temp_seurat <- subset(temp_seurat, nCount_RNA < 25000 & nFeature_RNA > 250 & percent.mt < 10)
 print(c(subject,dim(temp_seurat)))
 female_v2_stats_post_filter <- rbind(female_v2_stats_post_filter, temp_seurat@meta.data[,c("Sample", "nFeature_RNA", "nCount_RNA")])
 assign(paste0("Seurat_",subject), temp_seurat)
}

#Filter v3 female datasets by gene and UMI counts, recalculate gene and UMI stats
print("Filtering female v3 subjects")
for(subject in female_subs_v3$Sample) {
 temp_seurat <- get(paste0("Seurat_", subject))
 temp_seurat <- subset(temp_seurat, nCount_RNA < 120000 & nFeature_RNA > 350 & percent.mt < 10)
 print(c(subject,dim(temp_seurat)))
 female_v3_stats_post_filter <- rbind(female_v3_stats_post_filter, temp_seurat@meta.data[,c("Sample", "nFeature_RNA", "nCount_RNA")])
 assign(paste0("Seurat_",subject), temp_seurat)
}

#summarize gene and UMIs counts per dataset chemistry combo after filtering 
print("Summarizing post filtering stats for males")
summary(male_stats_post_filter)
print("Summarizing post filtering stats for v2 females")
summary(female_v2_stats_post_filter)
print("Summarizing post filtering stats for v3 females")
summary(female_v3_stats_post_filter)

save.image(file = "/home/malosree/scratch/Finalized_merging_normalization/1_Finalized_objs_preprocess_per_sub.RData")

print("Print warnings summary.")
summary(warnings())


print("Print session info.")
sessionInfo()

