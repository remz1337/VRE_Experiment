#!/bin/sh

EXPERIMENT_ID=$1

#PROXY_JOB_ID=$(sbatch --parsable sbatch_proxy.sh)
#PREPARE_JOB_ID=$(sbatch --parsable --dependency=after:$PROXY_JOB_ID+1 sbatch_prepare.sh $EXPERIMENT_ID)

PREPARE_JOB_ID=$(sbatch --parsable sbatch_prepare.sh $EXPERIMENT_ID)

#GENERATE_JOB_ID=$(sbatch --parsable --dependency=afterok:$PREPARE_JOB_ID sbatch_generate.sh $EXPERIMENT_ID)

echo "Launched experiment $EXPERIMENT_ID with slurm job $PREPARE_JOB_ID"
