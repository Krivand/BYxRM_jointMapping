#!/bin/bash
# can be called within project folder using ./mergeCSV.sh 
# don't forget to allow permissions using chmod +x mergeCSV.sh

# Define project path
PROJECTPATH=/home/albertf/vandy147/p05
IMED_PATH=${PROJECTPATH}/imed

# SPECIFY PARAMETERS
BASE_NAME="mergeCSV"
WALL_TIME="24:00:00"
MEMORY="4G"
MAIL_TYPE="ALL"
MAIL_USER="vandy147@umn.edu"

# Array of job parameters
declare -a JOBS=(
  "inputPath=${IMED_PATH}/joint_hotspot_QTLs outputPath=${IMED_PATH} rName=jhQTLs cleanUp=FALSE"
  "inputPath=${IMED_PATH}/unfiltered_assocs outputPath=${IMED_PATH} rName=unfilteredAssocs cleanUp=FALSE"
)

# Unset R_HOME to avoid conflicts
unset R_HOME

# Load the R module (make sure to specify the version you want)
module load R/4.4.0-openblas-rocky8

# Submit jobs for each job configuration
for JOB in "${JOBS[@]}"; do
  # Extract job-specific parameters
  eval $JOB

  sbatch --job-name=${rName} \
         --output=${PROJECTPATH}/logs/${rName}.out \
         --error=${PROJECTPATH}/logs/${rName}.err \
         --time=$WALL_TIME \
         --ntasks=1 \
         --cpus-per-task=1 \
         --mem=$MEMORY \
         --mail-type=$MAIL_TYPE \
         --mail-user=$MAIL_USER \
         --wrap="Rscript ${PROJECTPATH}/mergeCSV.R $inputPath $outputPath $rName $cleanUp"
done