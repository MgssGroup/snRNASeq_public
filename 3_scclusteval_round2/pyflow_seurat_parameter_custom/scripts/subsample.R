library(Seurat)
library(tidyverse)
library(harmony)

## see https://bitbucket.org/snakemake/snakemake/issues/917/enable-stdout-and-stderr-redirection
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

seurat_obj<- readRDS(snakemake@input[[1]])
k<- snakemake@wildcards[["k"]]
resolution<- snakemake@wildcards[["resolution"]]
pc.use<- snakemake@wildcards[["pc"]]
run_id<- snakemake@wildcards[["run_id"]]

PreprocessSubsetData_pars<- snakemake@params[["PreprocessSubsetData_pars"]]

RandomSubsetData<- function(object, rate, random.subset.seed = NULL, ...){
        ncells<- nrow(object@meta.data)
        ncells.subsample<- round(ncells * rate)

        set.seed(random.subset.seed)

        selected.cells<- sample(colnames(object), ncells.subsample)
        object<- subset(object, cells =  selected.cells,
                            ...)
        return(object)
}

subset_seurat_obj<- RandomSubsetData(seurat_obj, rate = snakemake@params[["rate"]])
original_ident<- Idents(subset_seurat_obj)

PreprocessSubsetData_custom <- function (object, num.pc = 20, pc.use = NULL,
    workers = 2, score.thresh = 1e-05, sig.pc.thresh = 0.05,
    n.start = 100, nn.eps = 0, resolution = 0.8, k.param = 30,
    ...)
{
    if (!is.null(pc.use)) {
        if (pc.use > num.pc) {
            stop("Specify the maximum pc.use number as less than or equal to the total num.pc")
        }
    }
    meta.data.colnames <- object@meta.data %>% colnames()
    vars.to.regress <- c("percent.mt", "nCount_RNA")
    vars.to.regress <- vars.to.regress[vars.to.regress %in% meta.data.colnames]
    object <- ScaleData(object, vars.to.regress = vars.to.regress, model.use = "poisson", use.umi = TRUE)
    object <- RunPCA(object = object, features = VariableFeatures(object = object),
        npcs = num.pc)
    pc.use.meta <- rep(pc.use, length(colnames(object)))
    names(pc.use.meta) <- colnames(object)
    object <- AddMetaData(object = object, metadata = pc.use.meta,
        col.name = "pc.use")
    object <- RunHarmony(object, reduction.save = "Harmony_Batch_Sample_Chemistry", 
			group.by.vars = c("Batch", "Sample", "Chemistry"))
    object <- UpdateSeuratObject(object)
    object <- FindNeighbors(object, reduction = "Harmony_Batch_Sample_Chemistry", dims = 1:pc.use, k.param = k.param,
        nn.eps = nn.eps, verbose = FALSE, force.recalc = TRUE)
    object <- FindClusters(object = object, n.start = n.start,
        resolution = resolution, verbose = FALSE, algorithm = 2)
    return(object)
}

## after reprocessing, the ident slot will be updated with the new cluster id
command<- paste("PreprocessSubsetData_custom", "(", "subset_seurat_obj,", "k.param=", k, ",", "pc.use=", pc.use, ",",
                                   "resolution=", resolution, ",", PreprocessSubsetData_pars, ")")
subset_seurat_obj<- eval(parse(text=command))

res<- tibble::tibble(pc = pc.use, resolution = resolution, k_param = k, original_ident = list(original_ident),
    recluster_ident = list(Idents(subset_seurat_obj)), round = run_id)


outfile<- paste0("subsample/subsample_", "k_", k, "_resolution_", resolution, "_PC_", pc.use, "_round_", run_id, ".rds")
saveRDS(res, file = outfile)

## make sure it is not empty file
info<- file.info(outfile)
if (info$size == 0) {
    quit(status = 1)
}

print(summary(warnings()))
print(sessionInfo())
