#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --gres=gpu:k80:1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=30G
#SBATCH --time=04:00:00
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/post_snakemake.err
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/post_snakemake.out

Rscript /home/malosree/projects/def-gturecki/malosree/scclusteval_round2/pyflow_seurat_parameter_custom/post_snakemake.R
