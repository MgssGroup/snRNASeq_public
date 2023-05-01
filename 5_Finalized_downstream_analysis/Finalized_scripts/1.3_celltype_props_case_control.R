library(Seurat)
library(tidyverse)
library(ggplot2)
library(rstatix)
library(boot)
library(parallel)

#Load data
print("Load data")
harmonized_all <- readRDS("/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/1_enhanced_harmonized_object.Rds")

cluster_cats <- c("Cluster", "Broad")

R <- 10000
parallel <- "no"
ncpus <- 1
#set.seed(21)

plots_list <- list()

for(sex in list(c("Male", "Female"), c("Female"), c("Male"))) {
set.seed(21)

harmonized_object <- subset(harmonized_all, Sex %in% sex)

#Get subject proportions in each cluster
print("Get subject proportions in each cluster")

for (cluster_cat in cluster_cats) {

subject_proportions <- table(harmonized_object$Sample, t(harmonized_object@meta.data[cluster_cat]))
if("Male" %in% sex) {
	subject_proportions["M24",] <- subject_proportions["M24",]+subject_proportions["M24_2",]
	subject_proportions <- subject_proportions[setdiff(rownames(subject_proportions), "M24_2"),]
}

#just a check
temp <- t(apply(subject_proportions, 1, function(x){
        x/sum(x)
}))

subject_proportions <- subject_proportions/rowSums(subject_proportions)

subject_proportions %>% as_tibble(.name_repair = "minimal") %>% setNames(c("Subject", "Cluster", "Proportion")) %>%
                pivot_wider(names_from = Cluster, values_from = Proportion) %>%
                mutate(Condition = map_chr(Subject, ~unique(harmonized_object$Condition[harmonized_object$Sample == .x]))) -> subject_proportions

#just a check
sum(temp == subject_proportions) 

#Perform wilcoxon tests for subject cell-type proprotions between cases and controls
#Similar to Habib et al., 2020 (Regev, Schwartz): https://doi.org/10.1038/s41593-020-0624-8

print("Perform wilcoxon tests for subject cell-type proprotions between cases and controls")
wilcoxon_proportions <- function(cluster_proportions, conditions){
	myData <- data.frame(Proportion = cluster_proportions, Condition = as.factor(conditions))
	res_test <- wilcox_test(formula = Proportion ~ Condition, data = myData, alternative = "two.sided", paired = FALSE, detailed = TRUE)
}

lapply(select(subject_proportions, where(is.numeric)), function(x) {
	wilcoxon_proportions(x, subject_proportions$Condition)
}) %>% bind_rows(.id = "Cluster")-> results_wilcoxon

#Perform bootstrapped wilcoxon tests for subject cell-type proportions between cases and controls
#Generate a bootstrapped p-value based on whether the estimate is different from 0 (greater or lesser)
#The estimate represents the median difference between values from the two conditions, as per the wilcox_test documentation
print("Perform bootstrapped wilcoxon tests for subject cell-type proportions between cases and controls")

w_test <- function(d, i){
  d <- d[i,]
  res <- wilcox_test(formula = Proportion ~ Condition, data = d, paired = FALSE, alternative = "two.sided", detailed = TRUE )
  return(res$estimate)	
}

lapply(select_if(subject_proportions, is.numeric), function(x){
        myData <- data.frame(Proportion = x, Condition = subject_proportions$Condition)                
	booted_res <- boot(myData, w_test, R=R, parallel = parallel, ncpus = ncpus)
	p1 <- 2 * sum(booted_res$t[,1] > 0) / R
	p2 <- 2 * sum(booted_res$t[,1] < 0) / R
  print(p1)
  print(p2)
  print(booted_res$t0) 
	p = ifelse(p1<1,p1,p2)
	tibble(estimate_booted = booted_res$t0, booted_p = p)
}) %>% bind_rows(.id = "Cluster") -> results_wilcoxon_booted


#Save the wilcoxon booted and unbooted results

results_wilcoxon <- full_join(results_wilcoxon, results_wilcoxon_booted, by = "Cluster")

results_wilcoxon %>% mutate(p_adjust = p.adjust(p, method = "BH")) %>% mutate(booted_p_adjust = p.adjust(booted_p, method = "BH")) -> results_wilcoxon

sex_label <- paste0(sex, collapse = "")

print(sex_label)

print(paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/1.3_wilcoxon_results",
                cluster_cat, "_", sex_label,
                ".csv"))

write.csv(results_wilcoxon, file = paste0("/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/1.3_wilcoxon_results", 
		cluster_cat, "_", sex_label,
		".csv"))

#Make box plots for cell-type proportions between cases and controls

for(name in names(select_if(subject_proportions, is.numeric))) {
	print(name)
	myData <- data.frame(Proportion = subject_proportions[name][[1]], Condition = as.factor(subject_proportions$Condition), 
			Subject = subject_proportions$Subject)
	res_plot <- ggplot(data = myData, aes(x= Proportion, y = Condition, label = Subject, color = Condition))+geom_boxplot()+geom_text(size = 7)+coord_flip()+
			ggtitle(paste0(name,"_", sex_label))+
			theme_classic(base_size = 20)
	plots_list[[paste0(name, sex_label)]] <- res_plot

}

}
}

#Print boxplots out to file
pdf(file = "/home/malosree/projects/def-gturecki/malosree/Finalized_downstream_analysis/Finalized_outputs/1.3_boxplots_wilcoxon.pdf",
        onefile = TRUE)

        plots_list

dev.off()

print(summary(warnings()))
sessionInfo()
