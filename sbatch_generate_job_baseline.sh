#!/bin/sh
##SBATCH --job-name=SINGLETON
#SBATCH --time=00-06:00:00
#SBATCH --account=slurm-account
#SBATCH --mail-user=myemail@email.com
### Too many jobs, don't email for successful ones (array should notify only for entire array)
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
##SBATCH --mem-per-cpu=12G
#SBATCH --mem=24G
#Too many logs
##SBATCH --output=slurm/generate_baseline/slurm-%j.out
##SBATCH --error=slurm/generate_baseline/slurm-%j.err
#SBATCH --output=/dev/null
##SBATCH --error=/dev/null

echo "Slurm job id:$SLURM_JOB_ID"
#echo "Slurm array job id:$SLURM_ARRAY_JOB_ID "
#echo "Task id:$SLURM_ARRAY_TASK_ID"
echo "Running Baseline Java programs in parallel"

###First generate population
#Load required modules
module load viennarna
module load java
module load mariadb

# Unset to use more than 2Gb
unset JAVA_TOOL_OPTIONS

#Start the Mysql Proxy
ROOT_FOLDER=$PWD
#PROXYSCRIPT="$ROOT_FOLDER/proxy.sh"
#source $PROXYSCRIPT

EXPERIMENT_ID=$1

LOAD_VARIABLES="$ROOT_FOLDER/load_variables.sh"
source $LOAD_VARIABLES

#Tools
#ROOT_FOLDER=$PWD
#DATABASE="$ROOT_FOLDER/vre/vre.db3"
#CLEANUP="$ROOT_FOLDER/utils/RemoveTemporaryFiles.sh"
#MYSQL_EXEC="$ROOT_FOLDER/mysql_exec.sh"
#MYSQL_EXEC_JOB="$ROOT_FOLDER/mysql_exec_job.sh"

#Common filenames
#MYSQL_CREDENTIALS_FILENAME="mysql_creds.env"
#FINALCSV_FILENAME="vre.csv"
#BNK_OSC_BE_NAME="ECJ_BNK_BE"
#BNK_OSC_VRE_NAME="ECJ_BNK_VRE"

#Load MySQL credentials
#MYSQL_CREDENTIALS_FILE="$ROOT_FOLDER/$MYSQL_CREDENTIALS_FILENAME"
#MYSQL_CONN="mysql --defaults-extra-file=$MYSQL_CREDENTIALS_FILE -t -N"

#Get the related walk ids first
#WALK_IDS=$($MYSQL_EXEC "select DISTINCT(j.fk_walk_id) from baseline_job as j inner join walk as w on w.ID=j.fk_walk_id inner join sample as s on s.ID=w.fk_sample_id INNER join experiment as e on e.ID=s.fk_experiment_ID where j.already_generated=0 and e.ID=$EXPERIMENT_ID;")
WALK_IDS=$($MYSQL_EXEC "select DISTINCT(j.fk_walk_id) from baseline_job as j inner join walk as w on w.ID=j.fk_walk_id inner join sample as s on s.ID=w.fk_sample_id INNER join experiment as e on e.ID=s.fk_experiment_ID where j.completed=0 and e.ID=$EXPERIMENT_ID;")

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

#If no more jobs left, exit
if [ ${#WALK_IDS_STR} -gt 0 ]; then
	#Pick a job ID
	#JOB_ID=$($MYSQL_EXEC "BEGIN; SET @selected_job_id := 0; UPDATE baseline_job SET already_generated = 1, ID = (SELECT @selected_job_id := ID) WHERE already_generated = 0 and fk_walk_id in ($WALK_IDS_STR) LIMIT 1; SELECT @selected_job_id; COMMIT;")
	JOB_ID=$($MYSQL_EXEC "BEGIN; SET @selected_job_id := 0; UPDATE baseline_job SET started=1, started_time=NOW(), slurm_job=$SLURM_JOB_ID, ID = (SELECT @selected_job_id := ID) WHERE started = 0 and fk_walk_id in ($WALK_IDS_STR) LIMIT 1; SELECT @selected_job_id; COMMIT;")
	#Selected job id should not be 0
	
	echo "JOB_ID:$JOB_ID"

	#if JOB_ID == 0, exit
	while [ $JOB_ID -gt 0 ]; do
		#First, pull job details from DB
		#JOB_DETAILS=$($MYSQL_CONN --execute="select j.folder, j.java_args, j.population_file from baseline_job as j where j.ID=$JOB_ID;")
		JOB_DETAILS=$($MYSQL_EXEC_JOB "select j.folder, j.java_args, j.population_file from baseline_job as j where j.ID=$JOB_ID;")

		#split variables, starts at 2 with awk to skip the header delimiter (+------------+-------...|)
		JOB_FOLDER=$(echo $JOB_DETAILS | awk -F "|" '{gsub(/\s$/, "", $2);gsub(/^\s/, "", $2); print $2}')
		JAVA_ARGS=$(echo $JOB_DETAILS | awk -F "|" '{gsub(/\s$/, "", $3);gsub(/^\s/, "", $3); print $3}')
		POP_FILE=$(echo $JOB_DETAILS | awk -F "|" '{gsub(/\s$/, "", $4);gsub(/^\s/, "", $4); print $4}')

		#First, check if it has some leftovers from a previous job that would've crashed and delete those files
		rm $POP_FILE > /dev/null 2>&1
		
		SLURM_LOG_FOLDER="$ROOT_FOLDER/slurm/generate_baseline"
		CRITICAL_LOGS="$SLURM_LOG_FOLDER/CriticalLog-$SLURM_JOB_ID.txt"

		#Launch Java
		#cd $JOB_FOLDER
		TMP_FOLDER="$SLURM_TMPDIR/$JOB_ID/$SLURM_ARRAY_TASK_ID"
		mkdir -p "$TMP_FOLDER"
		#cd $SLURM_TMPDIR
		cd $TMP_FOLDER
		TMP_FINALCSV_FILE="$TMP_FOLDER/$VRE_POPULATION_FILENAME"
		#FINALCSV_FILE="$JOB_FOLDER/$FINALCSV_FILENAME"

		#java $JAVA_ARGS
		#Write pop to temp folder, then copy back
		#java $JAVA_ARGS -p stat.VRE-CSV=$TMP_FINALCSV_FILE > /dev/null 2>&1
		java -Xmx20g $JAVA_ARGS -p stat.VRE-CSV=$TMP_FINALCSV_FILE
		
		#Get the exit code. 0 means no error, else something happened
		JAVA_STATUS=$?
		#echo $JAVA_STATUS

		cp "$TMP_FINALCSV_FILE" "$POP_FILE"

		#Check that it worked by reading at least 2 lines (header + 1 gen)
		#Else fail the job to run again
		#NO_LINES=$(wc -l $TMP_FINALCSV_FILE | awk '{ print $1 }')
		CP_NO_LINES=$(wc -l $POP_FILE | awk '{ print $1 }')
		
		if [ $CP_NO_LINES -gt 1 ] && [ $JAVA_STATUS -eq 0 ]; then
			$MYSQL_EXEC "UPDATE baseline_job SET completed=1, completed_time=NOW() WHERE ID = $JOB_ID;"	   
		else
			echo "Java exit code:$JAVA_STATUS" > $CRITICAL_LOGS
			echo "Copied file is too small (or inconsistent individual length)." >> $CRITICAL_LOGS
			echo "CP_NO_LINES:$CP_NO_LINES" >> $CRITICAL_LOGS
			echo "EXPECTED_POP_SIZE:$EXPECTED_POP_SIZE" >> $CRITICAL_LOGS
			echo "TARGET_SIZE:$TARGET_SIZE" >> $CRITICAL_LOGS
			echo "LAST_SIZE:$LAST_SIZE" >> $CRITICAL_LOGS
			
			#cp "$TMP_JAVA_LOGS" "$JAVA_LOGS"
			#cp "$TMP_JAVA_ERR_LOGS" "$JAVA_ERR_LOGS"
		fi

		#source $CLEANUP

		cd $ROOT_FOLDER

		#JOB_ID=$($MYSQL_EXEC "BEGIN; SET @selected_job_id := 0; UPDATE baseline_job SET already_generated = 1, ID = (SELECT @selected_job_id := ID) WHERE already_generated = 0 and fk_walk_id in ($WALK_IDS_STR) LIMIT 1; SELECT @selected_job_id; COMMIT;")
		JOB_ID=$($MYSQL_EXEC "BEGIN; SET @selected_job_id := 0; UPDATE baseline_job SET started=1, started_time=NOW(), slurm_job=$SLURM_JOB_ID, ID = (SELECT @selected_job_id := ID) WHERE started = 0 and fk_walk_id in ($WALK_IDS_STR) LIMIT 1; SELECT @selected_job_id; COMMIT;")
	done

fi

#Cancel any remaining pending job once all VRE jobs have been processed
#but first check if there are pending jobs to avoid calling scancel too many times
PENDING_SLURM_JOBS=$(sq --jobs=$SLURM_ARRAY_JOB_ID -h -t PENDING -r | wc -l)

if [ $PENDING_SLURM_JOBS -gt 0 ]; then
	echo "There are $PENDING_SLURM_JOBS pending jobs to cancel."
	scancel --state=PENDING $SLURM_ARRAY_JOB_ID
fi
