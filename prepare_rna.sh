#!/bin/sh

#This script is called from prepare_walk.sh

#Need to find BE output and create initial population file in ECJ format
#+pass the pop.file argument to ECJ parameters

#Use BNK_SOURCE_EXP_ID, SLURM_ARRAY_TASK_ID (1 to MAX_SAMPLES inclusive), walk_it (0 to MAX_WALKS exclusive) to compute SQL Query limit


#ROOT_FOLDER=$1
#echo "ROOT_FOLDER:$ROOT_FOLDER"

#First, retrieve SAMPLE_ID using SLURM_ARRAY_TASK_ID
echo "First, retrieve SAMPLE_ID using SLURM_ARRAY_TASK_ID"
#echo "SQL:select tmp.ID from (select s.* from sample as s where s.fk_experiment_ID=$BNK_SOURCE_EXP_ID order by s.ID ASC limit $SLURM_ARRAY_TASK_ID) as tmp order by tmp.ID desc limit 1;"
SOURCE_SAMPLE_ID=$($MYSQL_EXEC "select tmp.ID from (select s.* from sample as s where s.fk_experiment_ID=$BNK_SOURCE_EXP_ID order by s.ID ASC limit $SLURM_ARRAY_TASK_ID) as tmp order by tmp.ID desc limit 1;")
echo "SOURCE_SAMPLE_ID:$SOURCE_SAMPLE_ID"

#Retrieve source walk id
walk_limit=walk_it
#+1 to walk_it since it starts at 0
((walk_limit++))
#echo "SQL:select tmp.ID from (select w.* from walk as w where w.fk_sample_id=$SOURCE_SAMPLE_ID order by w.ID ASC limit $walk_limit) as tmp order by tmp.ID DESC limit 1;"
SOURCE_WALK_ID=$($MYSQL_EXEC "select tmp.ID from (select w.* from walk as w where w.fk_sample_id=$SOURCE_SAMPLE_ID order by w.ID ASC limit $walk_limit) as tmp order by tmp.ID DESC limit 1;")
echo "SOURCE_WALK_ID:$SOURCE_WALK_ID"


##### Get path to VRE population.in file and copy in target walk folder
SOURCE_EXPERIMENT_FOLDER="$ROOT_FOLDER/experiments/exp_$BNK_SOURCE_EXP_ID"
SOURCE_SAMPLE_FOLDER="$SOURCE_EXPERIMENT_FOLDER/sample_$SOURCE_SAMPLE_ID"
SOURCE_WALK_FOLDER="$SOURCE_SAMPLE_FOLDER/walk_$SOURCE_WALK_ID"

#source
SOURCE_INITIAL_ECJ_POPULATION_FILE="$SOURCE_WALK_FOLDER/$INITIAL_ECJ_POPULATION_FILENAME"
#echo "SOURCE_INITIAL_ECJ_POPULATION_FILE:$SOURCE_INITIAL_ECJ_POPULATION_FILE"

#target
#INITIAL_ECJ_POPULATION_FILE="$WALK_FOLDER/$INITIAL_ECJ_POPULATION_FILENAME"
#echo "TARGET_INITIAL_ECJ_POPULATION_FILE:$INITIAL_ECJ_POPULATION_FILE"

#copy
cp "$SOURCE_INITIAL_ECJ_POPULATION_FILE" "$INITIAL_ECJ_POPULATION_FILE"