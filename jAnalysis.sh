#!/bin/bash
# Can be called within the project folder using ./jAnalysis.sh
# Don't forget to allow permissions using chmod +x jAnalysis.sh

# Define project path
PROJECTPATH=/home/albertf/vandy147/p05

# SPECIFY PARAMETERS
BASE_NAME="jAnalysis"
MAIL_TYPE="ALL"
MAIL_USER="vandy147@umn.edu"
CPUS_PER_TASK=1

# Define job types, memory requirements, and time limits
declare -A JOB_MEMORIES=( ["analyte"]="10G" ["medium"]="50G" )
declare -A JOB_TIMES=( ["analyte"]="6:00:00" ["medium"]="12:00:00" )
declare -a JOB_TYPES=("analyte" "medium")

# Unset R_HOME to avoid conflicts
unset R_HOME

# Load the R module
module load R/4.4.0-openblas-rocky8

# Submit jobs for each job type
for JOB_TYPE in "${JOB_TYPES[@]}"; do
  MEMORY=${JOB_MEMORIES[$JOB_TYPE]}
  TIME=${JOB_TIMES[$JOB_TYPE]}

  sbatch --job-name=${BASE_NAME}_${JOB_TYPE} \
         --output=${PROJECTPATH}/logs/${BASE_NAME}_${JOB_TYPE}.out \
         --error=${PROJECTPATH}/logs/${BASE_NAME}_${JOB_TYPE}.err \
         --time=$TIME \
         --ntasks=1 \
         --cpus-per-task=$CPUS_PER_TASK \
         --mem=$MEMORY \
         --mail-type=$MAIL_TYPE \
         --mail-user=$MAIL_USER \
         --wrap="Rscript ${PROJECTPATH}/${BASE_NAME}.R ${JOB_TYPE}"
done