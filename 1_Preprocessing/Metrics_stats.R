#load libraries
library(ggbiplot)
library(dplyr)
library(readr)
library(ggplot2)
library(reshape2)
library(tidyr)


#load metrics
metrics_all <- read_csv(file = "C:/Users/Debajyoti Saha/Desktop/Turecki lab/Reanalysis_May_2021_selected_files/Preprocessing/Male_female_metadata_combined.csv")

#check column specs
spec(metrics_all) %>% print
# cols(
#   Sample = col_character(),
#   OriginalSub = col_character(),
#   Condition = col_character(),
#   Bank = col_character(),
#   Batch = col_character(),
#   Chemistry = col_character(),
#   Sequencing = col_character(),
#   Age = col_double(),
#   PMI = col_double(),
#   pH = col_double(),
#   Race = col_character(),
#   Sex = col_character(),
#   EstimatedNumCells = col_double(),
#   MeanReadsCell = col_double(),
#   MedianGenesCell = col_double(),
#   NumReads = col_double(),
#   ValidBarcodes = col_double(),
#   SequencingSaturation = col_double(),
#   Q30Barcode = col_double(),
#   Q30RNARead = col_double(),
#   Q30UMI = col_double(),
#   ReadsMappedGenome = col_double(),
#   ReadsMappedConfGenome = col_double(),
#   ReadsMappedIntergenic = col_double(),
#   ReadsMappedIntronic = col_double(),
#   ReadsMappedExonic = col_double(),
#   ReadsMappedTranscriptome = col_double(),
#   ReadsMappedAntisense = col_double(),
#   FractionReadsinCells = col_double(),
#   TotalGenesDetected = col_double(),
#   MedianUMICountsCell = col_double()
# )

#Create new variable combining Sex and Condition info
metrics_all %>% unite("SexCondition", Sex, Condition, remove = FALSE) -> metrics_all

#Find averages of metrics per group for Sex, Condition, and SexCondition
# metrics %>% 
#   group_by(Condition, Sex) %>%
#   select(where(is.numeric)) %>% 
#   summarise_all(.funs = ~ mean(.x, na.rm = TRUE)) -> means_sex_condition
# metrics %>% 
#   group_by(Condition) %>%
#   select(where(is.numeric)) %>% 
#   summarise_all(.funs = ~ mean(.x, na.rm = TRUE)) -> means_condition
# metrics %>% 
#   group_by(Sex) %>%
#   select(where(is.numeric)) %>% 
#   summarise_all(.funs = ~ mean(.x, na.rm = TRUE)) -> means_sex

#check if distributions are normal (perhaps this needs to be done broken down by group)
metrics_all %>% 
  select(where(is.numeric)) %>% 
  summarise_all(.funs = ~shapiro.test(na.omit(.))$p.value) -> shapiro_p.vlaues

t(shapiro_p.vlaues) %>% print
# [,1]
# Age                      1.411794e-03
# PMI                      1.034249e-03
# pH                       3.735143e-01
# EstimatedNumCells        4.214870e-06
# MeanReadsCell            1.331527e-07
# MedianGenesCell          3.287507e-03
# NumReads                 2.804001e-08
# ValidBarcodes            1.911810e-12
# SequencingSaturation     3.421129e-07
# Q30Barcode               2.005025e-08
# Q30RNARead               5.017490e-07
# Q30UMI                   1.052196e-08
# ReadsMappedGenome        1.419065e-05
# ReadsMappedConfGenome    9.404798e-05
# ReadsMappedIntergenic    3.474332e-01
# ReadsMappedIntronic      5.222698e-01
# ReadsMappedExonic        1.120317e-02
# ReadsMappedTranscriptome 2.462458e-03
# ReadsMappedAntisense     1.235162e-03
# FractionReadsinCells     3.130883e-04
# TotalGenesDetected       1.304022e-02
# MedianUMICountsCell      7.205996e-08

#check if distributions are normal broken down by group
# metrics %>% 
#   group_by(Condition, Sex) %>%
#   select(where(is.numeric)) %>%
#   summarise_all(.funs = ~shapiro.test(na.omit(.))$p.value) -> grouped_shapiro_p.vlaues

#for making qqplots
# for(name in names(metrics)){
#   if(is.numeric(metrics[[name]])){
#     qqnorm(metrics[[name]], pch = 1, frame = FALSE)
#     qqline(metrics[[name]], col = "steelblue", lwd = 2, sub = name)
#   }
# }

#Run KW test for each metric for groups split by SexCondition
KW_SexCondition <- vector()

for(name in names(metrics_all)){
  metrics <- metrics_all
  if(name %in% c("Age", "pH", "PMI")) {
    metrics %>% filter(Sample != "M24_2") -> metrics
  }
  if(is.numeric(metrics[[name]])){
    #data_plot <- metrics[,c(name, Condition, Sex)]    
    KW_res <- kruskal.test(metrics[[name]] ~ metrics$SexCondition)
    print(name)
    print(KW_res)
    KW_SexCondition <- rbind(KW_SexCondition, c(name, KW_res$statistic, KW_res$p.value))
  }
}

colnames(KW_SexCondition) <- c("Metric", "Stat", "P.Value")
write.csv(KW_SexCondition, file = "C:/Users/Debajyoti Saha/Desktop/Turecki lab/Reanalysis_May_2021_selected_files/Preprocessing/KW_SexCondition.csv")

#Run KW test (equivalent to Wilcoxon test) for each metric split by Condition
KW_Condition <- vector()

for(name in names(metrics_all)){
  metrics <- metrics_all
  if(name %in% c("Age", "pH", "PMI")) {
    metrics %>% filter(Sample != "M24_2") -> metrics
  }
  if(is.numeric(metrics[[name]])){
    #data_plot <- metrics[,c(name, Condition, Sex)]    
    KW_res <- kruskal.test(metrics[[name]] ~ metrics$Condition)
    print(name)
    print(KW_res)
    KW_Condition <- rbind(KW_Condition, c(name, KW_res$statistic, KW_res$p.value))
  }
}

colnames(KW_Condition) <- c("Metric", "Stat", "P.Value")
write.csv(KW_Condition, file = "C:/Users/Debajyoti Saha/Desktop/Turecki lab/Reanalysis_May_2021_selected_files/Preprocessing/KW_Condition.csv")

#Replace missing value for pH with mean. 
metrics_all$pH %>% replace_na(mean(metrics_all$pH, na.rm = TRUE)) -> metrics_all$pH

#reshape dataframe for plotting
metrics_melted <- melt(metrics_all)

#Plot boxplots of each metric split by SexCondition
pdf(file = "C:/Users/Debajyoti Saha/Desktop/Turecki lab/Reanalysis_May_2021_selected_files/Preprocessing/boxplots_metrics.pdf", onefile = TRUE, height = 8, width = 10)

(ggplot(metrics_melted, aes(y = value, x = Condition, color = Sex)) +
  facet_wrap(~variable, scales = "free_y", ncol = 4) +
  geom_boxplot() +
  theme_classic()) %>% print

#with subject labels 
(ggplot(metrics_melted, aes(y = value, x = Condition, color = Sex)) +
  facet_wrap(~variable, scales = "free_y", ncol = 4) +
  geom_boxplot() +
  geom_text(aes(label = OriginalSub), size = 2, position = position_jitter(width = 0.1)) +
  theme_classic()) %>% print

dev.off()

#Select all numeric metrics for PCA, sun PCA
metrics_all %>% select(where(is.numeric)) %>% as.data.frame() -> for_pca
rownames(for_pca) <- paste(metrics_all$Sample)

res_prcomp <- prcomp(for_pca, scale. = TRUE)

#Save PC1 PC2 plots colored by different groupings 
pdf(file = "C:/Users/Debajyoti Saha/Desktop/Turecki lab/Reanalysis_May_2021_selected_files/Preprocessing/pca_metrics.pdf", onefile = TRUE)

ggbiplot(res_prcomp, ellipse=TRUE, groups=metrics$Condition, var.axes = FALSE) %>% print
ggbiplot(res_prcomp, ellipse=TRUE, groups=metrics$Chemistry, var.axes = FALSE) %>% print
ggbiplot(res_prcomp, ellipse=TRUE, groups=metrics$Sex, var.axes = FALSE) %>% print
ggbiplot(res_prcomp, ellipse=TRUE, groups=metrics$SexCondition, var.axes = FALSE) %>% print
ggbiplot(res_prcomp, ellipse=TRUE, groups=metrics$Batch, var.axes = FALSE) %>% print

dev.off()
