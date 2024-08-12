#!/bin/bash
# Can be called within the project folder using ./iStudyTraces.sh
# Don't forget to allow permissions using chmod +x iStudyTraces.sh

# Define project path
PROJECTPATH=/home/albertf/vandy147/p05

# SPECIFY PARAMETERS
BASE_NAME="iStudyTraces"
CPUS_PER_TASK=1
WALL_TIME="8:00:00"
MAIL_TYPE="ALL"
MAIL_USER="vandy147@umn.edu"

# Unset R_HOME to avoid conflicts
unset R_HOME

# Load the R module
module load R/4.4.0-openblas-rocky8

# Define job names and their corresponding memory requirements
declare -A JOB_MEMORIES=( ["Analyte"]="10G" ["Medium"]="10G" )
declare -a JOBS=("Analyte" "Medium")

# Loop through the jobs and submit each one
for JOB_SUFFIX in "${JOBS[@]}"; do
  JOB_NAME="${BASE_NAME}${JOB_SUFFIX}"
  MEMORY=${JOB_MEMORIES[$JOB_SUFFIX]}

  sbatch --job-name=${JOB_NAME} \
         --output=${PROJECTPATH}/logs/${JOB_NAME}.out \
         --error=${PROJECTPATH}/logs/${JOB_NAME}.err \
         --time=$WALL_TIME \
         --ntasks=1 \
         --cpus-per-task=$CPUS_PER_TASK \
         --mem=$MEMORY \
         --mail-type=$MAIL_TYPE \
         --mail-user=$MAIL_USER \
         --wrap="Rscript ${PROJECTPATH}/${JOB_NAME}.R"
done