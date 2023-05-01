#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH --account=rrg-gturecki
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=40G
#SBATCH --ntasks=8
#SBATCH --output="/home/malosree/scratch/Permutation_outputs/Slurm_files/Step2_"$1".out"
#SBATCH --error="/home/malosree/scratch/Permutation_outputs/Slurm_files/Step2_"$1".err"

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Finalized_scripts/2_male_broad.R $1
EOT
