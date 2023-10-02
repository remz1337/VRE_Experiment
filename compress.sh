#!/bin/sh

#Define some global constants
ROOT_FOLDER=$PWD
EXPERIMENT_ID=$1
#PREPARE_LOAD="$ROOT_FOLDER/prepare_load.sh"

#source $PREPARE_LOAD $EXPERIMENT_ID

#EXPERIMENT_FOLDER="$ROOT_FOLDER/experiments/exp_$EXPERIMENT_ID"
EXPERIMENT_FOLDER="exp_$EXPERIMENT_ID"
#mkdir -p "$EXPERIMENT_FOLDER"
COMPRESS_FOLDER="$ROOT_FOLDER/experiments"

cd "$COMPRESS_FOLDER"

#First create tarball, then compress
#tar -Ipigz --remove-files -cf "expr_$EXPERIMENT_ID.tar.gz" "$EXPERIMENT_FOLDER"
tar --exclude='*.png' --exclude='*.pdf' -Ipigz -cf "expr_$EXPERIMENT_ID.tar.gz" "$EXPERIMENT_FOLDER"

#Copy to persistent storage (project folder)
#First check if test DB or Prod
MYSQL_CREDENTIALS_FILENAME="mysql_creds.env"
MYSQL_CREDENTIALS_FILE="$ROOT_FOLDER/$MYSQL_CREDENTIALS_FILENAME"

database=$(awk -F "=" '/database/ {print $2}' $MYSQL_CREDENTIALS_FILE)

source_file="$COMPRESS_FOLDER/expr_$EXPERIMENT_ID.tar.gz"
target_file="/home/myuser/projects/slurm-account/myuser/VRE_Exp_BU/expr_$EXPERIMENT_ID.tar.gz"

if [[ $database == *"test"* ]]; then
    #move to test folder
	target_file="/home/myuser/projects/slurm-account/myuser/VRE_Exp_BU/test_exp/expr_$EXPERIMENT_ID.tar.gz"
fi

#Comment to keep in scratch only if the project folder is full
#cp "$source_file" "$target_file"

#Remove local files?
#rm "$source_file"