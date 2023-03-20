#!/bin/sh

#Define some global constants
ROOT_FOLDER=$PWD
EXPERIMENT_ID=$1

LOAD_VARIABLES="$ROOT_FOLDER/load_variables.sh"
source $LOAD_VARIABLES


#PREPARE_LOAD="$ROOT_FOLDER/prepare_load.sh"
#SBATCH_PREPARE_SAMPLE="$ROOT_FOLDER/sbatch_prepare_sample.sh"
#SBATCH_GENERATE="$ROOT_FOLDER/sbatch_generate.sh"
#PREPARE_SAMPLE="$ROOT_FOLDER/prepare_sample.sh"
#PREPARE_WALK="$ROOT_FOLDER/prepare_walk.sh"
#MYSQL_EXEC="$ROOT_FOLDER/mysql_exec.sh"

source $PREPARE_LOAD $EXPERIMENT_ID

$MYSQL_EXEC "update experiment set folder = \"$EXPERIMENT_FOLDER\" where ID = $EXPERIMENT_ID;"

#sample_it=0
#while [ "$sample_it" -lt "$MAX_SAMPLES" ]; do
#	#Start walk in background
#	source $PREPARE_SAMPLE $EXPERIMENT_ID
#	((sample_it++))
#done

#Build slurm array str
SAMPLE_JOB_ARRAY="1-$MAX_SAMPLES"
echo "SAMPLE_JOB_ARRAY:$SAMPLE_JOB_ARRAY"

#Limit the number of simultaneous jobs (because of MySQL connections)
#SAMPLE_JOB_ARRAY="$SAMPLE_JOB_ARRAY"

#Ensure same variable name is used by calling script
SAMPLE_JOB_ID=$(sbatch --parsable --array=$SAMPLE_JOB_ARRAY $SBATCH_PREPARE_SAMPLE $EXPERIMENT_ID)

GENERATE_JOB_ID=$(sbatch --parsable --dependency=afterok:$SAMPLE_JOB_ID $SBATCH_GENERATE $EXPERIMENT_ID)

cd $ROOT_FOLDER

echo Done!
