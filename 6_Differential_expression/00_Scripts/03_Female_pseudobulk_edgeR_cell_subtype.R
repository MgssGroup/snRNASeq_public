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
  "/lustre07/scratch/senthil/Corina_snRNAseq_in_MDD_R006421/Analysis_report/Figures_and_Tables_all/Female_dataset/01_pseudobulk/02_cell_subtype"

#setwd(mainDir)
loop_parameters <-
  as.data.frame(
    read_csv(
      "/lustre07/scratch/senthil/Corina_snRNAseq_in_MDD_R006421/Analysis_report/Figures_and_Tables_all/Female_dataset/Table_loop_parameters_Female.csv"
    )
  )
loop_parameters <-
  loop_parameters[(loop_parameters$ClusterID == "cell_type"), ]
loop_parameters  <-
  loop_parameters %>% mutate(File_No = row_number())
loop_parameters$ClusterID <- c("cell_subtype")
#cell_ID <- c("ExN10_L46")
cell_ID <- c("Mic")
#FEMALE##
MC_raw_clustered_Female <-
  readRDS(file.path(mainDir, "MC_raw_clustered_Female.rds"))
#        readRDS(file.path(mainDir,"MC_raw_clustered_Female.rds"))

for (i in 6:9) {
#for (i in 6:nrow(loop_parameters)) {
  # i = 10
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
                             "Chemistry")])
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
  # sample_id = "sample_id"
  # group_id = "group_id"
  # celltype_id = "cell_type"
  # covariates = "Batch"
  #
  
  #aggData####
  pb <- aggregateData(
    sce,
    assay = "counts",
    fun = "sum",
    by = c("cluster_id", "sample_id")
  )
  
  #pbDS####
  # run DS analysis c("edgeR", "DESeq2", "limma-trend", "limma-voom")
  try(res <-
        pbDS(
          pb,
          design = design,
          min_cells = min_cell_n,
          contrast = contrast,
          filter = "none",
          method = "edgeR"
        ))
  # png(
  #   file.path(
  #     outDir,
  #     paste0(
  #       run_name_dir,
  #       "/Fig_top-20_DS_genes_in_Mic_pbHeatmap.png"
  #     )
  #   ),
  #   width = 15,
  #   height = 10,
  #   units = 'in',
  #   res = 300
  # )
  # top-20 DS genes for single cluster
  # <- pbHeatmap(sce, res, k = cell_ID)
  #print(p2)
  #dev.off()
  if (i == 10) {
    saveRDS(object = res, file.path(outDir, paste0(run_name_dir, "/res.rds")))
    saveRDS(object = sce, file.path(outDir, paste0(run_name_dir, "/sce.rds")))
  } else{
    saveRDS(object = res, file.path(outDir, paste0(run_name_dir, "/res.rds")))
  }
  # # pre-computed expression frequencies & append
  # frq <-
  #   calcExprFreqs(sce, assay = "counts", th = 0)
  # # one assay per cluster
  # assayNames(frq)
  # 
  # # expression frequencies by
  # # sample & group; 1st cluster:
  # head(assay(frq))
  
  #resDS####
  #output_master# append CPMs & expression frequencies
  res <-
    resDS(sce,
          res,
          cpm = TRUE,
          frq = TRUE,
          bind = "row") # tidy format
  
  
  
  res <-  res %>%
    mutate(
      Case_sampleNoneZero = rowSums(across(
        c(
          "F1.frq",
          "F11.frq",
          "F12.frq",
          "F14.frq",
          "F15.frq",
          "F16.frq",
          "F17.frq",
          "F18.frq",
          "F19.frq",
          "F2.frq",
          "F20.frq",
          "F25.frq",
          "F27.frq",
          "F28.frq",
          "F3.frq",
          "F4.frq",
          "F5.frq",
          "F6.frq",
          "F8.frq",
          "F9.frq"
        ),
        ~ .x == 0
      )),
      Control_sampleNoneZero = rowSums(across(
        c(
          "F10.frq",
          "F13.frq",
          "F21.frq",
          "F22.frq",
          "F23.frq",
          "F24.frq",
          "F26.frq",
          "F29.frq",
          "F30.frq",
          "F31.frq",
          "F32.frq",
          "F33.frq",
          "F34.frq",
          "F35.frq",
          "F36.frq",
          "F37.frq",
          "F38.frq",
          "F7.frq"
        ),
        ~ .x == 0
      )),
      
    ) %>%
    mutate(Case_n9 = Case_sampleNoneZero > 8, Control_n9 = Control_sampleNoneZero > 8)
  res$Method <- run_name_dir
  
  #res1 <- res[res$Case_n9& res$Control_n9 == "TRUE" ,]
  #Method_name####
  write_csv(res, file = file.path(
    outDir,
    paste0(
      run_name_dir,
      "/Table_raw_masterfile_cpm_Freq_DS_stats.csv"
    )
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
  top3 <- top3[c(1, 2, 79:92)]
  write_csv(top3, file.path(
    outDir,
    paste0(
      run_name_dir,
      "/Table_Top3_DS_per_cluster_id_FDR_5Pct.csv"
    )
  ))
  
  #upset_plot####
  df <-  tbl_fil[c(1, 2, 79:92)]
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
