#!/bin/bash
#SBATCH --job-name=iLODs                # Job name
#SBATCH --output=iLODs.out	 	        # Standard output and error log
#SBATCH --error=iLODs.err 	 	        # Error log
#SBATCH --time=24:00:00                 # Time limit hrs:min:sec
#SBATCH --ntasks=1                      # Number of tasks (usually 1 for a single R script)
#SBATCH --cpus-per-task=11              # Number of CPU cores per task
#SBATCH --mem=48G                       # Memory per node
#SBATCH --mail-type=ALL	          	    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=vandy147@umn.edu	# Where to send mail

# Unset R_HOME to avoid conflicts
unset R_HOME

# Load the R module (make sure to specify the version you want)
module load R/4.4.0-openblas-rocky8

# Run the R script
Rscript ~/p05/iStudyLODs.R
