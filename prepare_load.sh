#!/bin/sh

#---------------- Variables

#Define some global constants
ROOT_FOLDER=$PWD
EXPERIMENT_ID=$1

LOAD_VARIABLES="$ROOT_FOLDER/load_variables.sh"
source $LOAD_VARIABLES

#------------------ Fetch details from DB

#Check experiment ID in DB and retrieve experiment type
EXPERIMENT_DETAILS=$($MYSQL_EXEC "select type,fk_problem_id,fk_algorithm_id,max_samples,max_walks,neighborhood,bnk_source_exp_id,target_phenotype from experiment where ID=$EXPERIMENT_ID limit 1;")

EXPERIMENT_TYPE=$(echo $EXPERIMENT_DETAILS | awk '{print $1}')
PROB_ID=$(echo $EXPERIMENT_DETAILS | awk '{print $2}')
ALGO_ID=$(echo $EXPERIMENT_DETAILS | awk '{print $3}')
MAX_SAMPLES=$(echo $EXPERIMENT_DETAILS | awk '{print $4}')
MAX_WALKS=$(echo $EXPERIMENT_DETAILS | awk '{print $5}')
NEIGHBORHOOD=$(echo $EXPERIMENT_DETAILS | awk '{print $6}')
BNK_SOURCE_EXP_ID=$(echo $EXPERIMENT_DETAILS | awk '{print $7}')
TARGET_PHENOTYPE=$(echo $EXPERIMENT_DETAILS | awk '{print $8}')

if [ "$EXPERIMENT_TYPE" != "vre" ] && [ "$EXPERIMENT_TYPE" != "be" ] && [ "$EXPERIMENT_TYPE" != "br" ]; then
	echo "Unable to locate experiment ID in the database or invalid experiment type."
	exit 1
elif [ -z "$MAX_SAMPLES" ] || [ $MAX_SAMPLES -le 0 ]; then
	echo "Number of samples must be greater than 0."
	exit 1
elif [ -z "$MAX_WALKS" ] || [ $MAX_WALKS -le 0 ]; then
	echo "Number of walks must be greater than 0."
	exit 1
elif [ -z "$NEIGHBORHOOD" ] || [ $NEIGHBORHOOD -le 0 ]; then
	if [ "$EXPERIMENT_TYPE" = "vre" ]; then
		echo "Number of neighborhoods must be greater than 0."
		exit 1	
	fi
fi

if [ "$EXPERIMENT_TYPE" == "be" ]; then
	BASELINE_KEEP_BEST=true
fi


#Retrieve Algorithm details
ALGORITHM_DETAILS=$($MYSQL_EXEC "select name,population_size,generations,other_parameters from algorithm where ID=$ALGO_ID limit 1;")
ALGORITHM_NAME=$(echo $ALGORITHM_DETAILS | awk '{print $1}')
POPULATION_SIZE=$(echo $ALGORITHM_DETAILS | awk '{print $2}')
GENERATIONS=$(echo $ALGORITHM_DETAILS | awk '{print $3}')
ALGO_OTHER_PARAMS=$(echo $ALGORITHM_DETAILS | awk '{print $4}')
#Parse the other parameters. Following this format: PARAM1=VALUE1;PARAM2=VALUE2;...
IFS=';' read -r -a ALGO_OTHER_PARAMS_ARRAY <<< "$ALGO_OTHER_PARAMS"
unset IFS
ECJ_JAR=$ROOT_FOLDER/$(echo ${ALGO_OTHER_PARAMS_ARRAY[0]} | cut -f2 -d\=)
ECJ_PARAM_FILE=$ROOT_FOLDER/$(echo ${ALGO_OTHER_PARAMS_ARRAY[1]} | cut -f2 -d\=)
RNG_SEED=$(echo ${ALGO_OTHER_PARAMS_ARRAY[2]} | cut -f2 -d\=)

#Retrieve Problem details
PROBLEM_DETAILS=$($MYSQL_EXEC "select name,genotype_length,rna_alphabet_size,bnk_gates,bnk_inputs,other_parameters from problem where ID=$PROB_ID limit 1;")
PROB_NAME=$(echo $PROBLEM_DETAILS | awk '{print $1}')
GENOTYPE_LENGTH=$(echo $PROBLEM_DETAILS | awk '{print $2}')
RNA_ALPHABET_SIZE=$(echo $PROBLEM_DETAILS | awk '{print $3}')
BNK_GATES=$(echo $PROBLEM_DETAILS | awk '{print $4}')
BNK_INPUTS=$(echo $PROBLEM_DETAILS | awk '{print $5}')
PROB_OTHER_PARAMS=$(echo $PROBLEM_DETAILS | awk '{print $6}')

#Setup command line arguments for ECJ. Using $ IN stat file to specify relative path
ECJ_BASE_ARGUMENTS="-jar $ECJ_JAR -file $ECJ_PARAM_FILE -p seed.0=$RNG_SEED -p pop.subpop.0.species.genome-size=$GENOTYPE_LENGTH"

EXPERIMENT_FOLDER="$ROOT_FOLDER/experiments/exp_$EXPERIMENT_ID"
mkdir -p "$EXPERIMENT_FOLDER"

TMP_EXPERIMENT_FOLDER="$SLURM_TMPDIR/experiments/exp_$EXPERIMENT_ID"
mkdir -p "$TMP_EXPERIMENT_FOLDER"