#!/bin/sh

###First generate population
#Load required modules
#module load mariadb

#ROOT_FOLDER=$PWD
ROOT_FOLDER=$(dirname $(readlink -f $0))

MAX_RETRIES=100

SQL_QUERY=$1
#SQL_QUERY="select type,fk_problem_id,fk_algorithm_id,max_samples,max_walks from experiment where ID=9 limit 1;"
MYSQL_CREDENTIALS_FILENAME="mysql_creds.env"

#Load MySQL credentials
MYSQL_CREDENTIALS_FILE="$ROOT_FOLDER/$MYSQL_CREDENTIALS_FILENAME"
MYSQL_CONN="mysql --defaults-extra-file=$MYSQL_CREDENTIALS_FILE -s"

#SQL_RESULTS=$(echo $SQL_QUERY | $MYSQL_CONN)
#SQL_RESULTS=$(echo $SQL_QUERY | $MYSQL_COPROC)

mysql_query_execution="false"
retry_cnt=0
while [ "$mysql_query_execution" == "false" ] &&  [ $retry_cnt -lt $MAX_RETRIES ]; do
	SQL_RESULTS=$(echo $SQL_QUERY | $MYSQL_CONN 2>&1)
	#save your output
	#echo "Mysql query ran successfully"
	
	if [[ $SQL_RESULTS == ERROR* ]] ; then
		#echo "Its an error:$SQL_RESULTS"
		#Failed, try again after 15 seconds		
		sleep 15
	else
		#echo "Yay it worked:$SQL_RESULTS"
		mysql_query_execution="true"
	fi

	((retry_cnt++))
done

if [ ${#SQL_RESULTS} -gt 0 ]; then
	echo "$SQL_RESULTS"
fi

#EXPERIMENT_DETAILS=$($MYSQL_CONN --execute="select type,fk_problem_id,fk_algorithm_id,max_samples,max_walks from experiment where ID=$EXPERIMENT_ID limit 1;")
#$MYSQL_CONN --execute="insert into sample (fk_experiment_id, random_id) values ($EXPERIMENT_ID, \"$SAMPLE_RANDOM_ID\");" 2> "$SAMPLE_RANDOM_ID_LOG"
