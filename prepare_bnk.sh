#!/bin/sh

#This script is called from prepare_walk.sh

#Need to find BE output and create initial population file in ECJ format
#+pass the pop.file argument to ECJ parameters

#Use BNK_SOURCE_EXP_ID, SLURM_ARRAY_TASK_ID (1 to MAX_SAMPLES inclusive), walk_it (0 to MAX_WALKS exclusive) to compute SQL Query limit

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

#echo "SQL:select j.population_file from baseline_job as j where j.fk_walk_id=$SOURCE_WALK_ID;"
BE_POP_FILE=$($MYSQL_EXEC "select j.population_file from baseline_job as j where j.fk_walk_id=$SOURCE_WALK_ID;")
echo "BE_POP_FILE:$BE_POP_FILE"



VRE_POPULATION_FILE="$WALK_FOLDER/$VRE_POPULATION_FILENAME"
INITIAL_ECJ_POPULATION_FILE="$WALK_FOLDER/$INITIAL_ECJ_POPULATION_FILENAME"
#ECJ_POP_FILE="$WALK_FOLDER/$INITIAL_ECJ_POPULATION_FILENAME"
TMP_INITIAL_ECJ_POPULATION_FILE="$TMP_WALK_FOLDER/$INITIAL_ECJ_POPULATION_FILENAME"
#TMP_ECJ_POP_FILE="$TMP_BNK_FOLDER/$INITIAL_ECJ_POPULATION_FILENAME"
#FINAL_POPULATION_FILE="$WALK_FOLDER/$VRE_POPULATION_FILENAME"
#JAVA_LOG_FILE="$WALK_FOLDER/$JAVA_LOG_FILENAME"
#VALID_PHENOTYPES_FILE="$WALK_FOLDER/$VALID_PHENOTYPES_FILENAME"
VALID_PHENOTYPES_FILE="$TMP_WALK_FOLDER/$VALID_PHENOTYPES_FILENAME"



#Copy file to TMP dir
TMP_BNK_FOLDER="$TMP_WALK_FOLDER/$SOURCE_WALK_ID"
TMP_BE_POP_FILE="$TMP_BNK_FOLDER/$VRE_POPULATION_FILENAME"


mkdir -p "$TMP_BNK_FOLDER"
#cd $TMP_BNK_FOLDER
cp "$BE_POP_FILE" "$TMP_BE_POP_FILE"

#Keep only first and last lines
FIRSTLINE=$(head -n 1 $TMP_BE_POP_FILE)
LASTLINE=$(tail -n 1 $TMP_BE_POP_FILE)

#echo "FIRSTLINE:$FIRSTLINE"
#echo "LASTLINE:$LASTLINE"

echo $FIRSTLINE > $TMP_BE_POP_FILE
echo $LASTLINE >> $TMP_BE_POP_FILE

#Convert to ECJ format
#echo "Convert to ECJ format"
$VRE_ECJ_BUILDER -i $TMP_BE_POP_FILE -n > $TMP_INITIAL_ECJ_POPULATION_FILE

#check
#filecheck=$(head -n 2 $TMP_BE_POP_FILE)
#echo "filecheck:$filecheck"

#copy to persistent storage
#cp "$TMP_INITIAL_ECJ_POPULATION_FILE" "$INITIAL_ECJ_POPULATION_FILE"
cp "$TMP_INITIAL_ECJ_POPULATION_FILE" "$INITIAL_ECJ_POPULATION_FILE"

#Add pop to ECJ
#ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p pop.file=$INITIAL_ECJ_POPULATION_FILE"

#echo "Completed prepare! Exiting..."
#exit 1