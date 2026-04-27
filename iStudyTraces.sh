#!/bin/bash

# Project path
PROJECTPATH=/PHShome/kv1026/p05

# JOB SETUP
BASE_NAME="iStudyTraces"
CPUS_PER_TASK=1
WALL_TIME="8:00:00"
MAIL_TYPE="ALL"
MAIL_USER="krisna_vandyke@mail.dfci.harvard.edu"

# Ensure logs exist
mkdir -p "${PROJECTPATH}/logs"

# Clean R env and load R
unset R_HOME
module load R/4.4.0

# Job-specific memory
declare -A JOB_MEMORIES=( ["Analyte"]="10G" ["Medium"]="10G" )
declare -a JOBS=("Analyte" "Medium")

# Submit both jobs
for JOB_SUFFIX in "${JOBS[@]}"; do
  JOB_NAME="${BASE_NAME}${JOB_SUFFIX}"
  MEMORY=${JOB_MEMORIES[$JOB_SUFFIX]}

  sbatch --job-name="${JOB_NAME}" \
         --output="${PROJECTPATH}/logs/${JOB_NAME}.out" \
         --error="${PROJECTPATH}/logs/${JOB_NAME}.err" \
         --time="${WALL_TIME}" \
         --ntasks=1 \
         --cpus-per-task="${CPUS_PER_TASK}" \
         --mem="${MEMORY}" \
         --mail-type="${MAIL_TYPE}" \
         --mail-user="${MAIL_USER}" \
         --wrap="Rscript ${PROJECTPATH}/${JOB_NAME}.R"
done