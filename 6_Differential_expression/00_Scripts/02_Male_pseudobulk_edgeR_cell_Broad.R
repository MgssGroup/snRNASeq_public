closeAllConnections()
rm(list = ls())

library(ComplexHeatmap)
library(cowplot)
library(dplyr)
library(future)
library(ggeasy)
library(ggplot2)
library(limma)
library(muscat)
library(muscatWrapper)
library(purrr)
library(readr)
library(scater)
library(Seurat)
library(SingleCellExperiment)
library(UpSetR)

#narval /lustre07/scratch/senthil/
#cpu_n = 40
#statice /Users/skailasam/Desktop/mnt/projects/

cpu_n = 40

options(future.globals.maxSize = 90000 * 1024 ^ 2)
plan("multiprocess", workers = cpu_n)
plan()
mainDir <-
  "/lustre07/scratch/senthil/Corina_snRNAseq_in_MDD_R006421/Analysis_report/Figures_and_Tables_all"
setwd(mainDir)
outDir <-
  "/lustre07/scratch/senthil/Corina_snRNAseq_in_MDD_R006421/Analysis_report/Figures_and_Tables_all/Male_dataset/01_pseudobulk/01_cell_Broad"

#setwd(mainDir)
loop_parameters <-
  as.data.frame(
    read_csv(
      "/lustre07/scratch/senthil/Corina_snRNAseq_in_MDD_R006421/Analysis_report/Figures_and_Tables_all/Male_dataset/Table_loop_parameters_Male.csv"
    )
  )
loop_parameters <-
  loop_parameters[(loop_parameters$ClusterID == "cell_type"),]
loop_parameters  <-
  loop_parameters %>% mutate(File_No = row_number())
loop_parameters$SeuratObject <- c("Male")
#loop_parameters$ClusterID <- c("cell_subtype")
cell_ID <- c("Ast2")
#cell_ID <- c("Ast")


MC_raw_clustered_Male <-
  readRDS(file.path(mainDir, "MC_raw_clustered_Male.rds"))
# temp.df <- read_csv("../../../../Additional_metadata_20220119.csv")
# sample_ID <- MC_raw_clustered_Male@meta.data$Sample
# RIN <- temp.df$RIN[match(sample_ID,temp.df$Sample)]
# MC_raw_clustered_Male <- AddMetaData(MC_raw_clustered_Male, RIN, col.name = "RIN")

MC_raw_clustered_Male <-
  subset(MC_raw_clustered_Male, subset = Sample != c("M24_2"))
#        readRDS(file.path(mainDir,"MC_raw_clustered_Female.rds"))
#MC_raw_clustered_Male[["Chemistry"]] <- MC_raw_clustered_Male@meta.data$Chemistry

#saveRDS(MC_raw_clustered_Male,"MC_raw_clustered_Male.rds")
for (i in 1:nrow(loop_parameters))  {
  #i = 10
  run_name_dir <-
    paste(
      loop_parameters$File_No[i],
      loop_parameters$SeuratObject[i],
      loop_parameters$MethodType[i],
      loop_parameters$Tool[i],
      loop_parameters$Covariates[i],
      loop_parameters$ClusterID[i],
      sep = "_"
    )
  
  #run_name####
  print(run_name_dir)
  #nthreads####
  
  #setwd(mainDir)
  if (file.exists(file.path(outDir, run_name_dir))) {
    setwd(file.path(outDir, run_name_dir))
  } else {
    dir.create(file.path(outDir, run_name_dir), recursive =  T)
    setwd(file.path(outDir, run_name_dir))
    
  }
  
  if (loop_parameters$SeuratObject[i] == "Male") {
    sce <- as.SingleCellExperiment(MC_raw_clustered_Male)
  } else{
    sce <- as.SingleCellExperiment(MC_raw_clustered_Female)
  }
  #min_cell_param
  if (loop_parameters$ClusterID[i] == "cell_type") {
    min_cell_n <- c(10)
  } else{
    min_cell_n <- c(5)
  }
  
  sce <- sce[rowSums(counts(sce) > 0) > 0, ]
  #colData(sce$cell_subtype) <- as.character(colData(sce$cell_subtype))
  sce <- prepSCE(
    sce,
    kid = paste0(loop_parameters$ClusterID[i]),
    sid = "Sample",
    gid = "Condition",
    drop = FALSE
  )
  
  #contrast####
  #design & contrast matrix
  temp <- colData(sce)
  rownames(temp) <- NULL
  anno <-
    unique.data.frame(temp[c("sample_id",
                             "group_id",
                             "Sex",
                             "Batch",
                             "Age",
                             "pH",
                             "PMI",
                             "RIN")])
  anno$pH[anno$sample_id == "F35"] <- 6.47
  #anno <- metadata(sce)$experiment_info
  #design <- model.matrix(~ 0 + group_id + factor(Sex), data = anno)
  
  if (loop_parameters$Covariates[i] == "None") {
    design <-
      eval(parse(text = paste(
        "model.matrix(~ 0+group_id", " ,  data = anno)", sep = ""
      )))
  } else {
    design <-
      eval(parse(
        text = paste(
          "model.matrix(~ 0+group_id",
          loop_parameters$CovariatesPlus[i],
          " ,  data = anno)",
          sep = ""
        )
      ))
  }
  
  dimnames(design) <-
    list(anno$sample_id, levels = make.names(colnames(design)))
  design
  # Making a contrast for each pairwise comparison.
  contrast <-
    makeContrasts(
      Case_vs_Control = (group_idCase - group_idControl),
      levels = make.names(colnames(design))
    )
  sample_id = "sample_id"
  group_id = "group_id"
  celltype_id = "cell_type"
  covariates = "Batch"
  
  
  #aggData####
  pb <- aggregateData(
    sce,
    assay = "counts",
    fun = "sum",
    by = c("cluster_id", "sample_id")
  )
  
  #pbDS####
  # run DS analysis c("edgeR", "DESeq2", "limma-trend", "limma-voom")
  min_cell_n
  try(res <-
        pbDS(
          pb,
          design = design,
          min_cells = min_cell_n,
          contrast = contrast,
          filter = "none",
          method = "edgeR"
        ))
  
  
  #resDS####
  #output_master# append CPMs & expression frequencies
  try(res <-
        resDS(sce,
              res,
              cpm = TRUE,
              frq = TRUE,
              bind = "row"))
  
  
  # #pre-computed expression frequencies & append
  # frq <-
  #   calcExprFreqs(sce, assay = "counts", th = 0)
  # # one assay per cluster
  # assayNames(frq)
  # 
  #  # expression frequencies by
  #  # sample & group; 1st cluster:
  #  frq.df<- as.data.frame(assay(frq))
  #  colnames(frq.df) <- paste0(colnames(frq.df),".frq")
  # frq.df<-  frq.df %>%
  #   mutate(Case_sampleNoneZero = rowSums(across(c("M1.frq", "M10.frq", "M11.frq", "M14.frq", "M17.frq", "M18.frq",
  #                                                 "M23.frq", "M26.frq", "M28.frq", "M30.frq", "M32.frq", "M33.frq",
  #                                                 "M34.frq", "M4.frq", "M5.frq", "M6.frq", "M8.frq"), ~ .x > 0)),
  #          Control_sampleNoneZero = rowSums(across(c("M12.frq", "M13.frq", "M15.frq", "M16.frq", "M19.frq", "M2.frq",
  #                                                    "M20.frq", "M21.frq", "M22.frq", "M24.frq", "M27.frq",
  #                                                    "M29.frq", "M3.frq", "M31.frq", "M7.frq", "M9.frq"), ~ .x > 0)), ) %>%
  #   mutate(Case_n9=Case_sampleNoneZero > 8, Control_n9=Control_sampleNoneZero > 8  )
  # temp.frq <- data.frame(
  #   gene = rep(rownames(frq), length(assays(frq))),
  #   cluster_id = rep(assayNames(frq), each = nrow(frq)),
  #   do.call("rbind", as.list(assays(frq))),row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
  # frq <- .tidy(frq, ei, append = ".frq")
  # res <- inner_join(frq, res, by = c("gene", "cluster_id"))
  
  res <-  res %>%
    mutate(Case_sampleNoneZero = rowSums(across(c("M1.frq", "M10.frq", "M11.frq", "M14.frq", "M17.frq", "M18.frq",
                                                  "M23.frq", "M26.frq", "M28.frq", "M30.frq", "M32.frq", "M33.frq",
                                                  "M34.frq", "M4.frq", "M5.frq", "M6.frq", "M8.frq"), ~ .x > 0)),
           Control_sampleNoneZero = rowSums(across(c("M12.frq", "M13.frq", "M15.frq", "M16.frq", "M19.frq", "M2.frq",
                                                     "M20.frq", "M21.frq", "M22.frq", "M24.frq", "M27.frq",
                                                     "M29.frq", "M3.frq", "M31.frq", "M7.frq", "M9.frq"), ~ .x > 0)), ) %>%
    mutate(Case_n9=Case_sampleNoneZero > 8, Control_n9=Control_sampleNoneZero > 8  )
  
  # tidy format
  
  res$Method <- run_name_dir
  if (i == 10) {
    saveRDS(object = res, file.path(outDir, paste0(run_name_dir, "/res.rds")))
    saveRDS(object = sce, file.path(outDir, paste0(run_name_dir, "/sce.rds")))
  }else {
    
    saveRDS(object = res, file.path(outDir, paste0(run_name_dir, "/res.rds")))
  }
  #res1 <- res[res$Case_n9& res$Control_n9 == "TRUE" ,]
  #Method_name####
  write_csv(res, file = file.path(
    outDir,
    paste0(run_name_dir,
           "/Table_raw_masterfile_cpm_Freq_DS_stats.csv")
  ))
  #res<- read_csv("10_Female_Pseudobulk_edgeR_Covariate_Batch_Age_PMI_cell_type/Table_raw_masterfile_cpm_Freq_DS_stats.csv") %>% as.data.frame()
  
  
  # filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
  tbl_fil <-
    res %>% dplyr::filter(p_adj.loc < 0.05, abs(logFC) > 1) %>% dplyr::arrange(p_adj.loc)
  #for_comparison_Table
  
  
  # nb. of DS genes & % of total by cluster
  n_de <- as.vector(table(tbl_fil$cluster_id))
  col_n <-
    as.vector(names(table(tbl_fil$cluster_id)))
  p_de <-
    format(n_de / nrow(sce) * 100, digits = 3)
  Pct_DS <-
    cbind.data.frame("cell_ID" = col_n,
                     "n_DS" = n_de,
                     "Pct_DS" = p_de)
  write_csv(Pct_DS, file.path(
    outDir,
    paste0(
      run_name_dir,
      "/Table_DS_count_Proportion_per_cluster_id.csv"
    )
  ))
  
  
  # view top 2 hits in each cluster
  top3 <-
    tbl_fil %>% group_by(cluster_id) %>% arrange(p_adj.loc) %>% slice_head(n = 3)
  top3 <- top3[c(1, 2, 36:43)]
  write_csv(top3, file.path(
    outDir,
    paste0(run_name_dir,
           "/Table_Top3_DS_per_cluster_id_FDR_5Pct.csv")
  ))
  
  #upset_plot####
  df <-  tbl_fil[c(1, 2, 36:43)]
  write_csv(df, file.path(
    outDir,
    paste0(
      run_name_dir,
      "/Table_raw_masterfile_cpm_Freq_DS_stats_FDR_5Pct.csv"
    )
  ))
  
  de_gs_by_k <- split(df$gene, df$cluster_id)
  UpSet_fig <- upset(fromList(de_gs_by_k))
  
  
  png(
    file.path(
      outDir,
      paste0(run_name_dir, "/Fig_MDS_gene-overlap_UpSet_muscat.png")
    ),
    width = 15,
    height = 10,
    units = 'in',
    res = 300
  )
  
  print(UpSet_fig)
  #plot_grid(UpSet_fig, labels = c('A'))
  dev.off()
  
  
  # pull top-8 DS genes across all clusters
  top8 <- bind_rows(tbl_fil) %>%
    top_n(8, dplyr::desc(p_adj.loc)) %>%
    pull("gene")
  
  # for ea. gene in 'top8', plot t-SNE colored by its expression
  # wrapper to prettify reduced dimension plots
  .plot_dr <- function(sce, dr, col)
    plotReducedDim(sce, dimred = dr, colour_by = col) +
    guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_minimal() + theme(aspect.ratio = 1)
  
  png(
    file.path(
      outDir,
      paste0(run_name_dir, "/Fig_top-8_DS_genes_UMAP_colored.png")
    ),
    width = 15,
    height = 10,
    units = 'in',
    res = 300
  )
  ps <- lapply(top8, function(g)
    .plot_dr(sce, "UMAP_HARMONY_CHEMISTRY", g) +
      ggtitle(g) + theme(legend.position = "none"))
  
  
  # arrange plots 
  plot_grid(plotlist = ps,
            ncol = 4,
            align = "vh")
  dev.off()
  
  
}
