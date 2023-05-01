#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --error=/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Slurm_files/Step6.err
#SBATCH --output=/home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Slurm_files/Step6.out

Rscript /home/malosree/projects/def-gturecki/malosree/Finalized_cluster_QC_annotation/Finalized_scripts/6_cluster_naming.R
