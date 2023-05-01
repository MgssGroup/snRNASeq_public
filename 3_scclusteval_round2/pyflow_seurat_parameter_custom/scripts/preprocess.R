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

PreprocessFullData_custom <- function (object, num.pc = 20, pc.use = NULL,
    n.start = 100, nn.eps = 0, resolution = 0.8, k.param = 30,
    ...)
{
    if (!is.null(pc.use)) {
        if (pc.use > num.pc) {
            stop("Specify the maximum pc.use number as less than or equal to the total num.pc")
        }
    }
    pc.use.meta <- rep(pc.use, length(colnames(object)))
    names(pc.use.meta) <- colnames(object)
    object <- AddMetaData(object = object, metadata = pc.use.meta,
        col.name = "pc.use")
    object <- FindNeighbors(object, reduction = "Harmony_Batch_Sample_Chemistry", dims = 1:pc.use, k.param = k.param,
        nn.eps = nn.eps, verbose = FALSE, force.recalc = TRUE)
    object <- FindClusters(object = object, n.start = n.start,
        resolution = resolution, verbose = FALSE, algorithm = 2)
    return(object)
}

PreprocessSubsetData_pars<- snakemake@params[["PreprocessSubsetData_pars"]]
## this is not subsetted data, but the PreprocessSubsetData function can be used as well for any seurat object
seurat_obj<- eval(parse(text=paste("PreprocessFullData_custom", "(", "seurat_obj,", "k.param=", k, ",", "pc.use=", pc.use, ",",
                                   "resolution=", resolution, ",", PreprocessSubsetData_pars, ")")))
saveRDS(seurat_obj, file = paste0("full_sample_preprocess/full_sample_", "k_", k, "_resolution_", resolution, "_PC_", pc.use, ".rds"))

print(summary(warnings()))
print(sessionInfo())
