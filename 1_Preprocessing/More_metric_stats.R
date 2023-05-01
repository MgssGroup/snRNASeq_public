#load libraries
library(ggbiplot)
library(dplyr)
library(readr)
library(ggplot2)
library(reshape2)
library(tidyr)


#load metrics
metrics_all <- read_csv(file = "C:/Users/Debajyoti Saha/Desktop/Turecki lab/Reanalysis_May_2021_selected_files/Preprocessing/Male_female_metadata_combined.csv")

for(sex in c("Female", "Male")) {
  metrics_sex <- metrics_all %>% filter(Sex == sex)
  #check if distributions are normal (perhaps this needs to be done broken down by group)
  metrics_sex %>% 
    select(where(is.numeric)) %>% 
    summarise_all(.funs = ~shapiro.test(na.omit(.))$p.value) -> shapiro_p.vlaues
  print(t(shapiro_p.vlaues))
  #Run KW test for each metric for groups split by SexCondition
  
  KW_Condition <- vector()
  
  for(name in names(metrics_sex)){
    metrics <- metrics_sex
    if(name %in% c("Age", "pH", "PMI")) {
      metrics_sex %>% filter(Sample != "M24_2") -> metrics
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
  write.csv(KW_Condition, file = paste0("C:/Users/Debajyoti Saha//Desktop/Turecki lab/Reanalysis_May_2021_selected_files/Preprocessing/KW_", sex,".csv"))
  
}

#Shapiro female
# [,1]
# Age                      8.830498e-02
# PMI                      3.826144e-02
# pH                       5.176118e-01
# EstimatedNumCells        9.587838e-05
# MeanReadsCell            7.758039e-06
# MedianGenesCell          4.531861e-03
# NumReads                 3.456436e-09
# ValidBarcodes            3.498614e-07
# SequencingSaturation     1.232084e-02
# Q30Barcode               2.784723e-11
# Q30RNARead               2.892269e-04
# Q30UMI                   9.692356e-12
# ReadsMappedGenome        3.916245e-06
# ReadsMappedConfGenome    4.012774e-03
# ReadsMappedIntergenic    3.986652e-01
# ReadsMappedIntronic      1.197406e-01
# ReadsMappedExonic        3.818284e-01
# ReadsMappedTranscriptome 1.758387e-03
# ReadsMappedAntisense     2.421418e-03
# FractionReadsinCells     1.008550e-02
# TotalGenesDetected       4.462277e-02
# MedianUMICountsCell      3.555841e-05


#Shapiro male
# [,1]
# Age                      3.264891e-03
# PMI                      1.338960e-03
# pH                       2.859298e-01
# EstimatedNumCells        3.283761e-01
# MeanReadsCell            2.268604e-03
# MedianGenesCell          9.591379e-01
# NumReads                 7.338954e-01
# ValidBarcodes            1.085999e-09
# SequencingSaturation     7.225253e-01
# Q30Barcode               4.620802e-03
# Q30RNARead               1.750497e-01
# Q30UMI                   5.761899e-03
# ReadsMappedGenome        3.013320e-03
# ReadsMappedConfGenome    3.773045e-03
# ReadsMappedIntergenic    6.710148e-01
# ReadsMappedIntronic      2.654351e-01
# ReadsMappedExonic        5.708332e-02
# ReadsMappedTranscriptome 2.039114e-01
# ReadsMappedAntisense     7.745646e-01
# FractionReadsinCells     2.514559e-01
# TotalGenesDetected       4.134115e-05
# MedianUMICountsCell      9.983511e-01