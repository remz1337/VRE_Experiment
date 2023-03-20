#!/bin/sh
##SBATCH --job-name=SINGLETON
#SBATCH --time=00-06:00:00
#SBATCH --account=slurm-account
#SBATCH --mail-user=myemail@email.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
##SBATCH --mem=60G
#SBATCH --output=slurm/compress/slurm-%j.out
##SBATCH --error=slurm/compress/slurm-%j.err

EXPERIMENT_ID=$1
echo "Experiment id:$EXPERIMENT_ID"

echo "Slurm job id:$SLURM_JOB_ID"
#echo "Task id:$SLURM_ARRAY_TASK_ID"
echo "Compression job"

#ROOT_FOLDER=$PWD

#Execute program
./compress.sh $EXPERIMENT_ID
