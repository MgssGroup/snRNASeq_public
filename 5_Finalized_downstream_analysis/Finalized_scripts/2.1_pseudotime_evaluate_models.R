library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(tidyverse)
library(patchwork)
library(ggrastr)

print("Load data")
load(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/2_pseudotime.RData")

object_names <- c("male", "female")

#Replace missing pH value for a female subject with the average for the female dataset
average_female_pH <- as_tibble(colData(female)[c("Sample", "pH")]) %>% group_by(Sample) %>% summarise(pH = unique(pH)) %>%
		ungroup() %>% summarise(pH = mean(pH, na.rm = TRUE)) %>% unlist()
print(average_female_pH)
colData(female)$pH[is.na(colData(female)$pH)] <- rep(average_female_pH, sum(is.na(colData(female)$pH)))

#object_names <- c("female")

#object_names <- c("male")

BPPARAM <- BiocParallel::bpparam()

BPPARAM$workers <- 8

genes <- c("SOX10", "PDGFRA", "MYT1", "OLIG1", "OLIG2", "TCF7L2", "MBP", "PLP1", "GAPDH")

#Run fitGAM on range of ks and plot smoothers for marker genes of OL lineage
print("Run fitGAM on range of ks and plot smoothers for marker genes of OL lineage")
master_smoother_plots_list <- list()
#for(k in 3:20){
for(k in 3:10){
        smoother_plots_list <- list()
        for(object_name in object_names) {
                print(paste(object_name, k))
                sce <- get(object_name)
		print(paste("Fitting GAMS:", object_name, k))
                U_model <- model.matrix(~0+colData(sce)$Age+colData(sce)$PMI+colData(sce)$pH+colData(sce)$Batch)
                set.seed(20)
                sce <- fitGAM(counts = as.matrix(counts(sce)),
                           pseudotime = slingPseudotime(colData(sce)$slingshot),
                           cellWeights = slingCurveWeights(colData(sce)$slingshot),
                           conditions = factor(colData(sce)$Condition),
                           U = U_model,
                           nknots = k,
                           genes = genes,
                           parallel = TRUE,
                           BPPARAM = BPPARAM)
                print("Plotting smoothers")
                smoothers_individual <- list()
                for(gene in genes) {
                        print("Here")
                        smoothers_individual[[gene]] <- rasterize(plotSmoothers(sce, counts(sce), gene)+ggtitle(gene))
                }
                print("Here next")
		#Mistakenly changed and changed back on March 18, 2022
                patch <- Reduce("+", smoothers_individual) + plot_layout(ncol = 3, guides = "collect")
                smoother_plots_list[[object_name]] <- patch
        }
        master_smoother_plots_list[[as.character(k)]] <- smoother_plots_list
}

print("Saving smoother plots")
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Pseudotime/2.1_OLGene_smoothers.pdf", onefile = TRUE,
        height = 8, width = 12)
lapply(names(master_smoother_plots_list), function(x){
        master_smoother_plots_list[[x]]
})
dev.off()

print("Rerun evaluateK with 800 genes for male and female datasets separately.")
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/Pseudotime/2.1_evaluateK.pdf", 
				onefile = TRUE)
for(object_name in object_names) {
	for(seed in c(20, 99)) {
	print(paste(seed, object_name))
	sce <- get(object_name)
	sce <- sce[rowSums(counts(sce) > 0) >= 10, ]
	U_model <- model.matrix(~0+colData(sce)$Age+colData(sce)$PMI+colData(sce)$pH+colData(sce)$Batch)
	set.seed(seed)
	icMAT <- evaluateK(counts = as.matrix(counts(sce)),
                   pseudotime = slingPseudotime(colData(sce)$slingshot),
                   cellWeights = slingCurveWeights(colData(sce)$slingshot),
                   conditions = factor(colData(sce)$Condition),
                   U = U_model,
		   nGenes = 800,
                   #k = 3:20, 
		   k = 3:10,
		   parallel = TRUE,
		   BPPARAM = BPPARAM)
	}
 }
dev.off()

print(summary(warnings()))
sessionInfo()
