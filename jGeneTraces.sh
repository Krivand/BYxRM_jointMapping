#!/bin/bash
# Can be called within project folder using ./jGeneTraces.sh
# Don't forget to allow permissions using chmod +x jGeneTraces.sh

# Define project path
PROJECTPATH=/home/albertf/vandy147/p05

# SPECIFY PARAMETERS
BASE_NAME="jGeneTraces"
MAIL_TYPE="ALL"
MAIL_USER="vandy147@umn.edu"
TIME="8:00:00"
CPUS_PER_TASK=1

# Define job types and memory requirements
declare -A JOB_MEMORIES=( ["Analyte"]="10G" ["Medium"]="100G" )
declare -a JOB_TYPES=("Analyte" "Medium")

# Unset R_HOME to avoid conflicts
unset R_HOME

# Load the R module
module load R/4.4.0-openblas-rocky8

# Submit jobs for each job type
for JOB_TYPE in "${JOB_TYPES[@]}"; do
  MEMORY=${JOB_MEMORIES[$JOB_TYPE]}

  sbatch --job-name=${BASE_NAME}_${JOB_TYPE} \
         --output=${PROJECTPATH}/logs/${BASE_NAME}_${JOB_TYPE}.out \
         --error=${PROJECTPATH}/logs/${BASE_NAME}_${JOB_TYPE}.err \
         --time=$TIME \
         --ntasks=1 \
         --cpus-per-task=$CPUS_PER_TASK \
         --mem=$MEMORY \
         --mail-type=$MAIL_TYPE \
         --mail-user=$MAIL_USER \
         --wrap="Rscript ${PROJECTPATH}/${BASE_NAME}${JOB_TYPE}.R"
done