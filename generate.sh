#!/bin/sh

#---------------- Variables

#Define some global constants
ROOT_FOLDER=$PWD
EXPERIMENT_ID=$1

LOAD_VARIABLES="$ROOT_FOLDER/load_variables.sh"
source $LOAD_VARIABLES

#INVERSERNA_NAME="IMAP"
#SYSTEM_INTEGER_LIMIT=2147483647
#MAX_GENERATIONS_WITH_BUFFER=2045222520
#it's 10000, but using 9000 just in case
#MAX_SLURM_ARRAY_SIZE=9000

#Common filenames
#INITIAL_VRE_POPULATION_FILENAME="rnavre.csv"
#FINAL_POPULATION_FILENAME="vre_pop.csv"
#VALID_PHENOTYPES_FILENAME="valid_phenotypes.txt"
#INITIAL_ECJ_POPULATION_FILENAME="population.in"
#SQLITE_ERROR_FILENAME="sqlite_errors.err"
#MATH_HELPER_FILENAME="vre_math_helper.txt"

#Tools
#DATABASE="$ROOT_FOLDER/vre/vre.db3"
#VRE="$ROOT_FOLDER/vre/VRE"
#VRE_ECJ_BUILDER="$ROOT_FOLDER/vre/VRE_ECJ_Builder"
#VRNA_PHENOTYPE_GENERATOR="$ROOT_FOLDER/utils/RNA_PhenotypeGenerator.sh"
#VRE_MATH_HELPER="$ROOT_FOLDER/utils/VRE_MathHelper"
#VRE_PHENOTYPE_GENERATOR="$ROOT_FOLDER/utils/VRE_PhenotypeGenerator"
#SBATCH_JAVA_VRE="$ROOT_FOLDER/sbatch_generate_job_vre.sh"
#SBATCH_JAVA_BASELINE="$ROOT_FOLDER/sbatch_generate_job_baseline.sh"
#MYSQL_EXEC="$ROOT_FOLDER/mysql_exec.sh"
#PREPARE_LOAD="$ROOT_FOLDER/prepare_load.sh"

#ECJ Variables
#RNG_SEED=time
#BASELINE_KEEP_BEST=false

#RNAInverse
#INVERSE_REPEAT=50

source $PREPARE_LOAD $EXPERIMENT_ID

#------------------ Fetch details from DB

#Check experiment ID in DB and retrieve experiment type
#EXPERIMENT_DETAILS=$($MYSQL_EXEC "select type,fk_problem_id,fk_algorithm_id,max_samples,max_walks from experiment where ID=$EXPERIMENT_ID limit 1;")
#EXPERIMENT_TYPE=$(echo $EXPERIMENT_DETAILS | awk '{print $1}')
#PROB_ID=$(echo $EXPERIMENT_DETAILS | awk '{print $2}')
#ALGO_ID=$(echo $EXPERIMENT_DETAILS | awk '{print $3}')
#MAX_SAMPLES=$(echo $EXPERIMENT_DETAILS | awk '{print $4}')
#MAX_WALKS=$(echo $EXPERIMENT_DETAILS | awk '{print $5}')

#if [ "$EXPERIMENT_TYPE" != "vre" ] && [ "$EXPERIMENT_TYPE" != "be" ] && [ "$EXPERIMENT_TYPE" != "br" ]; then
#	echo "Unable to locate experiment ID in the database or invalid experiment type."
#	exit 1
#elif [ -z "$MAX_SAMPLES" ] || [ $MAX_SAMPLES -le 0 ]; then
#	echo "Number of samples must be greater than 0."
#	exit 1
#elif [ -z "$MAX_WALKS" ] || [ $MAX_WALKS -le 0 ]; then
#	echo "Number of walks must be greater than 0."
#	exit 1
#fi

if [ "$EXPERIMENT_TYPE" == "vre" ]; then
	JOB_TABLE="vre_neighborhood_job"
	SBATCH_JAVA=$SBATCH_JAVA_VRE
else
	JOB_TABLE="baseline_job"
	SBATCH_JAVA=$SBATCH_JAVA_BASELINE
fi

#Just count the number of jobs to run and call job array
PENDING_JOB_COUNT=$($MYSQL_EXEC "select count(j.ID) as job_count from $JOB_TABLE as j inner join walk as w on w.ID=j.fk_walk_id inner join sample as s on s.ID=w.fk_sample_id inner join experiment as e on e.ID=s.fk_experiment_ID where e.ID=$EXPERIMENT_ID and j.started=0;")

#Build slurm array str
JOB_ID_ARRAY="1-$PENDING_JOB_COUNT"
#JOB_ID_ARRAY="1-$MAX_SAMPLES"

#Limit the number of simultaneous jobs (because of MySQL connections AND CC user QoS limit)
if [ $PENDING_JOB_COUNT -gt $MAX_SLURM_ARRAY_SIZE ]; then
	JOB_ID_ARRAY="1-$MAX_SLURM_ARRAY_SIZE"
fi

JOB_ID_ARRAY="$JOB_ID_ARRAY"

#Dynamically adjust time needs, based on population size and number of generation
JOB_TIME="00-02:00:00"
EXP_SCALE=$((POPULATION_SIZE * GENERATIONS))

echo "GENERATIONS:$GENERATIONS"
echo "POPULATION_SIZE:$POPULATION_SIZE"
echo "EXP_SCALE:$EXP_SCALE"
echo "JOB_TIME:$JOB_TIME"

#Only need to scale for VRE
#if [ "$EXPERIMENT_TYPE" == "vre" ]; then
#	if [ $EXP_SCALE -ge 1000000 ]; then 
#		JOB_TIME="01-06:00:00"
#	elif [ $EXP_SCALE -ge 100000 ]; then 
#		JOB_TIME="00-12:00:00"
#	elif [ $EXP_SCALE -ge 10000 ]; then 
#		JOB_TIME="00-04:00:00"
#	else
#		JOB_TIME="00-08:00:00"
#	fi
#fi

if [ $EXP_SCALE -ge 1000000 ]; then 
	JOB_TIME="01-12:00:00"
elif [ $EXP_SCALE -ge 100000 ]; then 
	JOB_TIME="01-00:00:00"
elif [ $EXP_SCALE -ge 10000 ]; then 
	JOB_TIME="00-06:00:00"
else
	JOB_TIME="00-08:00:00"
fi


echo "JOB_TIME:$JOB_TIME"

#Ensure same variable name is used by calling script
GENERATE_JOB_ID=$(sbatch --parsable --time=$JOB_TIME --array=$JOB_ID_ARRAY $SBATCH_JAVA $EXPERIMENT_ID)