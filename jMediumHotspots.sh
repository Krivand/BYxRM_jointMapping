#!/bin/bash
# Can be called within the project folder using ./jMediumHotspots.sh 
# Don't forget to allow permissions using chmod +x jMediumHotspots.sh

# Define project path
PROJECTPATH=/home/albertf/vandy147/p05

# SPECIFY PARAMETERS
NUM_JOBS=42
CPUS_PER_TASK=64
BASE_NAME="jMediumHotspots"
MEMORY="120G"
MAIL_TYPE="ALL"
MAIL_USER="vandy147@umn.edu"

# Calculate the number of elements per job
TOTAL_ELEMENTS=84

# Adjust NUM_JOBS if it exceeds TOTAL_ELEMENTS
if [ $NUM_JOBS -gt $TOTAL_ELEMENTS ]; then
  NUM_JOBS=$TOTAL_ELEMENTS
fi

ELEMENTS_PER_JOB=$((TOTAL_ELEMENTS / NUM_JOBS))
REMAINDER=$((TOTAL_ELEMENTS % NUM_JOBS))

# Function to calculate wall time dynamically based on elements
calculate_wall_time() {
  local num_elements=$1
  local hours=$((12 * num_elements))
  printf "%02d:00:00" "$hours"
}

# Unset R_HOME to avoid conflicts
unset R_HOME

# Load the R module (make sure to specify the version you want)
module load R/4.4.0-openblas-rocky8

START=1

for ((i=0; i<NUM_JOBS; i++)); do
  if [ $i -lt $REMAINDER ]; then
    END=$((START + ELEMENTS_PER_JOB))
  else
    END=$((START + ELEMENTS_PER_JOB - 1))
  fi

  # Calculate wall time for this job
  NUM_ELEMENTS=$((END - START + 1))
  WALL_TIME=$(calculate_wall_time "$NUM_ELEMENTS")

  # Submit the job only if START is less than or equal to TOTAL_ELEMENTS
  if [ $START -le $TOTAL_ELEMENTS ]; then
    echo "Submitting job ${BASE_NAME}_${i} for range ${START} to ${END}" # Debug message
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
  fi

  START=$((END + 1))
done