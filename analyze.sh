#!/bin/sh

#---------------- Variables
ROOT_FOLDER=$PWD
EXPERIMENT_ID=$1

LOAD_VARIABLES="$ROOT_FOLDER/load_variables.sh"
source $LOAD_VARIABLES

#Common filenames
#MYSQL_CREDENTIALS_FILENAME="mysql_creds.env"

#Executables
#VRE="$ROOT_FOLDER/vre/VRE"
#MYSQL_EXEC="$ROOT_FOLDER/mysql_exec.sh"
#PREPARE_LOAD="$ROOT_FOLDER/prepare_load.sh"

#Load MySQL credentials
#MYSQL_CREDENTIALS_FILE="$ROOT_FOLDER/$MYSQL_CREDENTIALS_FILENAME"

#EXPERIMENT_FOLDER="$ROOT_FOLDER/experiments/exp_$EXPERIMENT_ID"

source $PREPARE_LOAD $EXPERIMENT_ID
#Dynamically adjust memory needs, based on population size and number of generation
PLOTTER_MEM="4G"
EXP_SCALE=$((POPULATION_SIZE * GENERATIONS))
if [ $EXP_SCALE -ge 1000000 ]; then 
	PLOTTER_MEM="160G"
elif [ $EXP_SCALE -ge 100000 ]; then 
	PLOTTER_MEM="80G"
elif [ $EXP_SCALE -ge 10000 ]; then 
	PLOTTER_MEM="40G"
else
	PLOTTER_MEM="16G"
fi


if [ "$EXPERIMENT_TYPE" == "vre" ]; then
	JOB_TABLE="vre_neighborhood_job"
	#SBATCH_JAVA=$SBATCH_JAVA_VRE
else
	JOB_TABLE="baseline_job"
	#SBATCH_JAVA=$SBATCH_JAVA_BASELINE
fi

#First, check if all jobs completed successfully
REMAINING_JOBS=$($MYSQL_EXEC "SELECT count(j.ID) as remaining_jobs FROM $JOB_TABLE as j inner join walk as w on w.ID=j.fk_walk_id inner join sample as s on s.ID=w.fk_sample_id inner join experiment as e on e.ID=s.fk_experiment_ID where e.ID=$EXPERIMENT_ID and j.completed=0;")
echo "There are $REMAINING_JOBS remaining jobs to process."

if [ $REMAINING_JOBS -eq 0 ]; then

	#VRE_ARGUMENTS="--experiment=$EXPERIMENT_ID --credentials=$MYSQL_CREDENTIALS_FILE -f"
	VRE_ARGUMENTS="--experiment=$EXPERIMENT_ID --credentials=$MYSQL_CREDENTIALS_FILE"

	#Run the VRE Analysis
	$VRE $VRE_ARGUMENTS

	#Then create the plots
	#sbatch --dependency=afterok:$ANALYZE_JOB_ID sbatch_plotter.sh $EXPERIMENT_ID
	#PLOTTER_JOB_ID=$(sbatch --parsable --dependency=afterok:$ANALYZE_JOB_ID sbatch_plotter.sh $EXPERIMENT_ID)
	PLOTTER_JOB_ID=$(sbatch --mem=$PLOTTER_MEM --parsable sbatch_plotter.sh $EXPERIMENT_ID)

	#Finally, compress and delete files
	COMPRESS_JOB_ID=$(sbatch --parsable --dependency=afterok:$PLOTTER_JOB_ID sbatch_compress.sh $EXPERIMENT_ID)
else
	#Missed some jobs, run again
	echo "Launching generate script again to complete missed jobs"

	#first reset started states where not completed
	WALK_IDS=$($MYSQL_EXEC "select DISTINCT(j.fk_walk_id) from $JOB_TABLE as j inner join walk as w on w.ID=j.fk_walk_id inner join sample as s on s.ID=w.fk_sample_id INNER join experiment as e on e.ID=s.fk_experiment_ID where j.completed=0 and e.ID=$EXPERIMENT_ID;")
	
	echo "Retrieved walks:$WALK_IDS"

	#Build string of walk ids
	WALK_IDS_STR=""
	walk_cnt=0
	while IFS= read -r WALK_ID; do
		#Append all the submitted job for next sbatch dependency
		if [ ${#WALK_ID} -gt 0 ]; then
			if [ $walk_cnt -eq 0 ]; then
				WALK_IDS_STR="$WALK_ID"
			else
				WALK_IDS_STR="$WALK_IDS_STR,$WALK_ID"
			fi
			((walk_cnt++))
	#	else
	#		echo "Empty walk ID:$WALK_ID, for experiment $EXPERIMENT_ID"
		fi
	done <<< "$WALK_IDS"

	$MYSQL_EXEC "UPDATE $JOB_TABLE SET started=0, started_time=NULL, slurm_job=NULL WHERE started=1 and completed=0 and fk_walk_id in ($WALK_IDS_STR);"
	
	echo "Reset job status complete."

	#echo "Temporarily disable automatic requeue. Go check the java log and understand how to capture failure earlier."
	#GENERATE_JOB_ID=$(sbatch --parsable sbatch_generate.sh $EXPERIMENT_ID)
	
	echo "Queued job generation again. Slurm ID:$GENERATE_JOB_ID"
fi
