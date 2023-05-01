#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH --account=def-naguibm
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --ntasks=8
#SBATCH --output="/home/malosree/projects/def-cnagy/malosree/Permutation_outputs/Slurm_files/Step4_"$1".out"
#SBATCH --error="/home/malosree/projects/def-cnagy/malosree/Permutation_outputs/Slurm_files/Step4_"$1".err"

module load r/4.1.2

Rscript /home/malosree/projects/def-gturecki/malosree/Permutation_analysis/Finalized_scripts/4_female_broad.R $1
EOT
