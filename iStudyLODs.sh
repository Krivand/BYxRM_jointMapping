#!/bin/bash
#SBATCH --job-name=iStudyLODs
#SBATCH --output=/PHShome/kv1026/p05/logs/iStudyLODs.out
#SBATCH --error=/PHShome/kv1026/p05/logs/iStudyLODs.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem=48G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=krisna_vandyke@mail.dfci.harvard.edu

# Clean environment
unset R_HOME

# Load R
module load R/4.4.0

# Run the R script
Rscript /PHShome/kv1026/p05/iStudyLODs.R
