library(scclusteval)
library(tidyverse)

subsample_idents_list <- readRDS(file =
        "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/subsample_idents_list.Rds")

fullsample_idents <- readRDS(file =
  "/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/gather_full_sample.rds")

fullsample_idents %>%
  mutate(cluster_num = purrr::map_dbl(original_ident_full, ~n_distinct(.x))) -> fullsample_idents

pdf(file = "/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_outputs/17_regen_ParameterScatters_70.pdf",
      height = 10, width = 12, onefile = TRUE)
subsample_idents_list_sub <- subsample_idents_list %>% mutate(stable_cluster = stable_cluster_70_0.7) %>% select(-res)
p1 <- ParameterSetScatterPlot(stable_clusters = subsample_idents_list_sub,
                        fullsample_idents = fullsample_idents,
                        x_var = "resolution",
                        y_var = "percentage",
                        facet_rows = "k_param",
                        facet_cols = "pc")
print(p1)
p1$data %>% select(resolution, k_param, pc, percentage) -> p1
write.csv(p1, file = 
	"/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_outputs/17_source_data_regen_ParameterScatters_70.csv")
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
        mutate(n_clusters = length(unlist(median_jaccards))) %>%
	#mutate(n_clusters = map(data, ~unique(n_distinct(.x$original_ident[[1]])))) %>%
	#old incorrect
	#mutate(n_clusters = map(data, ~n_distinct(.x$original_ident))) %>%
        unite(res, resolution, k_param, pc, remove = FALSE, sep = "_")

#create meds object to help with chooseR type plots
#slightly changed from before 
print("Creating meds object")
meds <- subsample_idents_list %>% select(n_clusters, boot, res, median_jaccards) %>% unnest_wider(boot) 
	#%>% 
	#mutate(n_clusters = unlist(n_clusters))

#Minimum thershold for stable clusters based on bootstrapping
#copied from chooseR
print("Calculating threshold")
# Find thresholds
threshold <- max(meds$low_med)
print(threshold)
choice <- as.character(
  meds %>%
  dplyr::filter(med >= threshold) %>%
  dplyr::arrange(n_clusters) %>%
  tail(n = 1) %>%
  dplyr::pull(res)
)

plot_data <- subsample_idents_list %>% select(res, median_jaccards) %>% unnest(median_jaccards)

# And plot!
print("Plotting chooseR type plot")
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_outputs/17_chooseR_bootstrapped_stability.pdf")
p1 <- ggplot(meds, aes(factor(res), med)) +
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
p1
p1$data %>% unnest(median_jaccards) -> p1_data
write.csv(p1_data, file = 
	"/home/malosree/projects/def-gturecki/malosree/Followup_analyses/Finalized_outputs/17_source_data_chooseR_bootstrapped_stability.csv")
dev.off()

print(summary(warnings()))
sessionInfo()
