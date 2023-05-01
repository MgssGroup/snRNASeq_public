#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH --account=def-naguibm
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=70G
#SBATCH --ntasks=8
#SBATCH --output="/home/malosree/projects/def-cnagy/malosree/Permutation_outputs/Slurm_files/Step1_"$1".out"
#SBATCH --error="/home/malosree/projects/def-cnagy/malosree/Permutation_outputs/Slurm_files/Step1_"$1".err"

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Finalized_scripts/1_male_subtype.R $1
EOT
