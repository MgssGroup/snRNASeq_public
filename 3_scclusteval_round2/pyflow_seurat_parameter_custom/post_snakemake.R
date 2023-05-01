library(scclusteval)
library(tidyverse)
library(patchwork)
library(Seurat)
library(dplyr)

print("Loading data")
#load clustering results for subsampled data and full data
subsample_idents <- readRDS(file = "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/gather_subsample.rds")

#prepare input for AssignStableCluster
subsample_idents_list<- subsample_idents %>%
  group_by(pc, resolution, k_param) %>%
  nest()

#free up some space 
rm(subsample_idents)

#find stable clusters
print("Finding stable clusters with stringent parameters")
subsample_idents_list <- subsample_idents_list %>%
  mutate(stable_cluster_80_0.8 = map(data, ~ AssignStableCluster(.x$original_ident,
                                                          .x$recluster_ident,
                                                          jaccard_cutoff = 0.8,
                                                          method = "jaccard_percent",
                                                          percent_cutoff = 0.8)))

print("Finding stable clusters with lenient parameters")
subsample_idents_list <- subsample_idents_list %>% 
	mutate(stable_cluster_70_0.7 = map(data, ~ AssignStableCluster(.x$original_ident,
                                                          .x$recluster_ident,
                                                          jaccard_cutoff = 0.7,
                                                          method = "jaccard_percent",
                                                          percent_cutoff = 0.7)))

print("Loading data")
#saveRDS(subsample_idents_list, file =
#      "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/subsample_idents_list.Rds")

fullsample_idents <- readRDS(file =
  "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/gather_full_sample.rds")

fullsample_idents %>%
  mutate(cluster_num = purrr::map_dbl(original_ident_full, ~n_distinct(.x))) -> fullsample_idents

print("Plotting")
pdf(file = "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/ParameterScatters.pdf", 
      height = 10, width = 12, onefile = TRUE)
subsample_idents_list <- subsample_idents_list %>% mutate(stable_cluster = stable_cluster_80_0.8)
ParameterSetScatterPlot(stable_clusters = subsample_idents_list,
                        fullsample_idents = fullsample_idents,
                        x_var = "resolution",
                        y_var = "number",
                        facet_rows = "k_param",
                        facet_cols = "pc")
ParameterSetScatterPlot(stable_clusters = subsample_idents_list,
                        fullsample_idents = fullsample_idents,
                        x_var = "resolution",
                        y_var = "percentage",
                        facet_rows = "k_param",
                        facet_cols = "pc")
dev.off()

pdf(file = "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/ParameterScatters_70.pdf",
      height = 10, width = 12, onefile = TRUE)
subsample_idents_list <- subsample_idents_list %>% mutate(stable_cluster = stable_cluster_70_0.7)
ParameterSetScatterPlot(stable_clusters = subsample_idents_list,
                        fullsample_idents = fullsample_idents,
                        x_var = "resolution",
                        y_var = "number",
                        facet_rows = "k_param",
                        facet_cols = "pc")
ParameterSetScatterPlot(stable_clusters = subsample_idents_list,
                        fullsample_idents = fullsample_idents,
                        x_var = "resolution",
                        y_var = "percentage",
                        facet_rows = "k_param",
                        facet_cols = "pc")
dev.off()

#boot median from chooseR

boot_median <- function(x, interval = 0.95, R = 25000, type = "bca") {
  # Define median to take data and indices for use with boot::
  med <- function(data, indices) {
    resample <- data[indices]
    return(median(resample))
  }

  # Calculate intervals
  boot_data <- boot::boot(data = x, statistic = med, R = R)
  boot_ci <- boot::boot.ci(boot_data, conf = interval, type = type)

  # Extract desired statistics
  ci <- list(
    low_med = boot_ci$bca[4],
    med = boot_ci$t0,
    high_med = boot_ci$bca[5]
  )
  return(ci)
}

#need to set seeds?

print("Calculating medians and bootstrapping")
#For each parameter set calculate median Jaccard indices per cluster
#For each parameter set bootstrap the median of the median Jaccard index per cluster
subsample_idents_list <- subsample_idents_list %>%
	mutate(stable_cluster = stable_cluster_80_0.8) %>%
        mutate(median_jaccards = map(stable_cluster, ~ unlist(lapply(.x$jaccardIndex, median)))) %>%
        mutate(boot = map(median_jaccards, boot_median)) %>%
	mutate(n_clusters = map(data, ~n_distinct(.x$original_ident))) %>%
	unite(res, resolution, k_param, pc, remove = FALSE, sep = "_")

#create meds object to help with chooseR type plots
print("Creating meds object")
meds <- subsample_idents_list %>% select(n_clusters, boot, res, median_jaccards) %>% unnest_wider(boot)

#Minimum thershold for stable clusters based on bootstrapping
#copied from chooseR
print("Calculating threshold")
# Find thresholds
threshold <- max(meds$low_med)
choice <- as.character(
  meds %>%
  dplyr::filter(med >= threshold) %>%
  dplyr::arrange(n_clusters) %>%
  tail(n = 1) %>%
  dplyr::pull(res)
)

plot_data <- subsample_idents_list %>% select(res, median_jaccards) %>% unnest(median_jaccards)

#copied from chooseR examples 1_seurat_pipeline.R and modified
# And plot!
print("Plotting chooseR type plot")
pdf(file = "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/chooseR_bootstrapped_stability.pdf")
ggplot(meds, aes(factor(res), med)) +
  geom_crossbar(
    aes(ymin = low_med, ymax = high_med),
    fill = "grey",
    size = 0.25
  ) +
  geom_jitter(
    data = plot_data,
    aes(factor(res), median_jaccards),
    size = 0.35,
    width = 0.15
  ) +
  geom_hline(aes(yintercept = threshold), colour = "blue") +
  geom_vline(aes(xintercept = choice), colour = "red") +
  scale_x_discrete("Parameter sets") +
  scale_y_continuous(
    "Jaccard Index",
    expand = c(0, 0),
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.25),
    oob = scales::squish
  ) +
  cowplot::theme_minimal_hgrid() +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7, angle = 90),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
  )
dev.off()

#count clusters that are stable based on the bootstrapped threshold and plot via scclusteval
print("Setting stable clusters by chooseR cut-off")
subsample_idents_list <- subsample_idents_list %>%
	mutate(stable_cluster_boot = map(data, ~ AssignStableCluster(.x$original_ident,
                                                          .x$recluster_ident,
                                                          jaccard_cutoff = threshold,
                                                          method = "jaccard_median")))

print("Plotting")
pdf(file = "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/ParameterScatters_boot.pdf",
      height = 10, width = 12, onefile = TRUE)
subsample_idents_list <- subsample_idents_list %>%
        mutate(stable_cluster = stable_cluster_boot)
#These plotting functions may not be working very well with the "stable_cluster_boot" ad-hoc column. They generate some warnings. 
ParameterSetScatterPlot(stable_clusters = subsample_idents_list,
                        fullsample_idents = fullsample_idents,
                        x_var = "resolution",
                        y_var = "number",
                        facet_rows = "k_param",
                        facet_cols = "pc")
ParameterSetScatterPlot(stable_clusters = subsample_idents_list,
                        fullsample_idents = fullsample_idents,
                        x_var = "resolution",
                        y_var = "percentage",
                        facet_rows = "k_param",
                        facet_cols = "pc")
dev.off()

subsample_idents_list %>% ungroup() %>% mutate(id = row_number()) %>%
   filter(pc == 70, resolution == 0.7, k_param == 30) -> temp
temp
temp$id
#[1] 23


pdf(file = "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/BestJaccard.pdf", 
	height = 10, width = 12, onefile = TRUE)
JaccardRainCloudPlot(subsample_idents_list$data[[8]]$original_ident,
                          subsample_idents_list$data[[8]]$recluster_ident) + 
        geom_hline(yintercept = c(0.7, 0.8, threshold), linetype = 2) +
        xlab("cluster id w/ k=30 res=0.7 pc=70") 
dev.off()

print("Saving data")
saveRDS(subsample_idents_list, file = "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/subsample_idents_list.Rds")

print(summary(warnings()))
sessionInfo()
