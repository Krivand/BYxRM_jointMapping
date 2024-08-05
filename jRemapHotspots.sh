#!/bin/bash
# can be called within project folder using ./jRemapHotspots.sh 
# don't forget to allow permissions using chmod +x jRemapHotspots.sh

# Define project path
PROJECTPATH=/home/albertf/vandy147/p05

# SPECIFY PARAMETERS
NUM_JOBS=8
BASE_NAME="jRemapHotspots"
WALL_TIME="24:00:00"
MEMORY="120G"
MAIL_TYPE="ALL"
MAIL_USER="vandy147@umn.edu"

# Calculate the number of elements per job
TOTAL_ELEMENTS=102
ELEMENTS_PER_JOB=$((TOTAL_ELEMENTS / NUM_JOBS))
REMAINDER=$((TOTAL_ELEMENTS % NUM_JOBS))

# Unset R_HOME to avoid conflicts
unset R_HOME

# Load the R module (make sure to specify the version you want)
module load R/4.4.0-openblas-rocky8

START=1

for ((i=0; i<NUM_JOBS; i++)); do
  END=$((START + ELEMENTS_PER_JOB - 1))
  if [ $i -lt $REMAINDER ]; then
    END=$((END + 1))
  fi

  CPUS_PER_TASK=$((END - START + 1))

  sbatch --job-name=${BASE_NAME}_${i} \
         --output=${PROJECTPATH}/logs/${BASE_NAME}_${i}.out \
         --error=${PROJECTPATH}/logs/${BASE_NAME}_${i}.err \
         --time=$WALL_TIME \
         --ntasks=1 \
         --cpus-per-task=$CPUS_PER_TASK \
         --mem=$MEMORY \
         --mail-type=$MAIL_TYPE \
         --mail-user=$MAIL_USER \
         --wrap="Rscript ${PROJECTPATH}/${BASE_NAME}.R ${START} ${END}"

  START=$((END + 1))
done