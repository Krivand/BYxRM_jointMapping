#!/bin/bash

# Define project path
PROJECTPATH=/home/albertf/vandy147/p05
BASE_NAME="jAnalysis"

#SBATCH --job-name=${BASE_NAME}                         # Job name
#SBATCH --output=${PROJECTPATH}/logs/${BASE_NAME}.out   # Standard output and error log
#SBATCH --error=${PROJECTPATH}/logs/${BASE_NAME}.err    # Error log
#SBATCH --time=8:00:00                                  # Time limit hrs:min:sec
#SBATCH --ntasks=1                                      # Number of tasks (usually 1 for a single R script)
#SBATCH --cpus-per-task=1                               # Number of CPU cores per task
#SBATCH --mem=8G                                        # Memory per node
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=vandy147@umn.edu                    # Where to send mail

# Unset R_HOME to avoid conflicts
unset R_HOME

# Load the R module (make sure to specify the version you want)
module load R/4.4.0-openblas-rocky8

# Run the R script
Rscript ${PROJECTPATH}/${BASE_NAME}.R