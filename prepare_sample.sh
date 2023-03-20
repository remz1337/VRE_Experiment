#!/bin/sh

#Define some global constants
ROOT_FOLDER=$PWD
EXPERIMENT_ID=$1
PREPARE_LOAD="$ROOT_FOLDER/prepare_load.sh"

source $PREPARE_LOAD $EXPERIMENT_ID

good_sample_random_id="false"
while [ "$good_sample_random_id" == "false" ]; do
	#Generate random id
	SAMPLE_RANDOM_ID=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1)
	
	SAMPLE_RANDOM_ID_LOG="$SLURM_TMPDIR/$SAMPLE_RANDOM_ID.log"

	rm $SAMPLE_RANDOM_ID_LOG > /dev/null 2>&1
	#Get a new Sample ID from SQLite DB
	$MYSQL_EXEC "insert into sample (fk_experiment_id, random_id) values ($EXPERIMENT_ID, \"$SAMPLE_RANDOM_ID\");" 2> "$SAMPLE_RANDOM_ID_LOG"

	if [[ -z $(grep '[^[:space:]]' $SAMPLE_RANDOM_ID_LOG) ]] ; then
		good_sample_random_id="true"
	fi
	
	#Cleanup
	rm $SAMPLE_RANDOM_ID_LOG > /dev/null 2>&1		
done

SAMPLE_ID=$($MYSQL_EXEC "select ID from sample where fk_experiment_id=$EXPERIMENT_ID and random_id=\"$SAMPLE_RANDOM_ID\";")

SAMPLE_FOLDER="$EXPERIMENT_FOLDER/sample_$SAMPLE_ID"
mkdir -p "$SAMPLE_FOLDER"

TMP_SAMPLE_FOLDER="$TMP_EXPERIMENT_FOLDER/sample_$SAMPLE_ID"
mkdir -p "$TMP_SAMPLE_FOLDER"

$MYSQL_EXEC "update sample set folder = \"$SAMPLE_FOLDER\" where ID = $SAMPLE_ID;"

echo "MAX_WALKS:$MAX_WALKS"

walk_it=0
while [ "$walk_it" -lt "$MAX_WALKS" ]; do

	echo "walk_it:$walk_it"

	cd $ROOT_FOLDER

	#Start walk in background
	source $PREPARE_WALK $EXPERIMENT_ID $SAMPLE_ID
	((walk_it++))
	
	echo "done walk!"
done
