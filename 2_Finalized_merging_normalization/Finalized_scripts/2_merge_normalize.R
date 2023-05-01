#load required libraries and seurat objects
library(Seurat)

load(file = "/home/malosree/scratch/Finalized_merging_normalization/1_Finalized_objs_preprocess_per_sub.RData")

object_list <- list()

#create and save object list, remove individual objects
for(object in ls(pattern = "^Seurat")){
 object_list[[object]] <- get(object)}

rm(list = ls(pattern = "^Seurat"))

#SCTransform each object in list, no covariates
print("Running SCTransform")
object_list <- lapply(X = object_list, FUN = function(x) {
    x <- SCTransform(x, method = "glmGamPoi") 
})

saveRDS(object_list, file = "/home/malosree/scratch/Finalized_merging_normalization/2_Finalized_preprocessed_objects_list.Rds")

#find variable features to be used for PCA later
var_features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000)

#make objects smaller, remove SCT assay and make DietSeurat
object_list <- lapply(X = object_list, FUN = function(x) {
   x@active.assay <- "RNA"
   x[["SCT"]] <- NULL
   x <- DietSeurat(x)
})

print("Merging objects")
#Merge objects 
objects_merged <- merge(object_list[[1]], y = object_list[2:length(object_list)], project = "combined", merge.data=FALSE)

rm(object_list)

#LogNormalize and assign the variable features. 
library(future)

plan("multicore", workers = 6)

options(future.globals.maxSize= 60*1000*1000*1000)

print("Normalizing")
objects_merged <- NormalizeData(objects_merged, normalization.method = "LogNormalize", scale.factor = 10000)

VariableFeatures(objects_merged) <- var_features

#Scale data regressing (linear) UMIs counts and percent.mito

print("Scaling")
objects_merged <- ScaleData(objects_merged, vars.to.regress = c("nCount_RNA", "percent.mt"), model.use = "poisson", use.umi = TRUE)

saveRDS(objects_merged, file = "/home/malosree/scratch/Finalized_merging_normalization/2_Finalized_merged_norm_obj.Rds")

print("Print warnings summary.")
summary(warnings())

print("Print session info.")
sessionInfo()

