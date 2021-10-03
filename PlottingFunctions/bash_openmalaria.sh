#!/bin/bash

#SBATCH --job-name=run_OpenMalaria
#SBATCH --account=pothin
#SBATCH --time=2:00:00
#SBATCH --qos=6hours
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=scicore_output.txt

############################################################
# BASH OPENMALARIA
#
# Submit OpenMalaria simulations as jobs to the cluster.
#
# Arguments:
#		  INPUT_DIR: Directory containing the scenario xml files
#		   DEST_DIR: Directory where OpenMalaria output files will be saved
#	 RESOURCE_DIR: Directory of OpenMalaria resources
#      LOG_FILE: Text log file that stores names of completed scenarios
#
# Created by M.Golumbeanu
# Adapted by A.J.Shattock
############################################################

# Load OpenMalaria module
module purge
ml OpenMalaria/38.0-goolf-1.7.20-Python-2.7.11

# Define input variables
INPUT_DIR=$1
DEST_DIR=$2
RESOURCE_DIR=$3
LOG_FILE=$4

# Change to folder with resource files
cd $RESOURCE_DIR

# Extract scenario xml files
SCENARIO_FILES=(${INPUT_DIR}*.xml)

# Select single scenario file based on task ID
ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)
SCENARIO_FILE=${SCENARIO_FILES[$ID]}

# Construct paths to xml file and output dirs
SCENARIO_NAME=$(basename "$SCENARIO_FILE" ".xml")
OUTPUT=$DEST_DIR$SCENARIO_NAME".txt"

# Check if the output file already exists
if [ -f "$OUTPUT1" ]; then
    echo "Scenario $SCENARIO_NAME has already been run - not re-running"
    
else 
    echo "Running simulation for $SCENARIO_NAME"

    # Run OpenMalaria for this scenario file
    openMalaria --scenario $SCENARIO_FILE --output $OUTPUT
fi

# Write scenario name to log file
echo "OpenMalaria simulation $SCENARIO_NAME finished" >> $LOG_FILE

