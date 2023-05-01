#!/bin/bash
#SBATCH --account=def-gturecki
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=50g
#SBATCH --time=24:00:00

module load mugqic/cellranger/5.0.1

cd /home/malosree/scratch/realign_May2020/male_counts

FASTQS=$1
NAME=$2

cellranger count --id=$NAME --fastqs=$FASTQS --transcriptome=/home/malosree/projects/def-gturecki/refdata-gex-GRCh38-2020-A/ --include-introns
