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
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#Too many logs
#SBATCH --output=slurm/prepare_sample/slurm-%j.out
##SBATCH --error=slurm/prepare_sample/slurm-%j.err
##SBATCH --output=/dev/null
##SBATCH --error=/dev/null

EXPERIMENT_ID=$1
echo "Experiment id:$EXPERIMENT_ID"

echo "Slurm job id:$SLURM_JOB_ID"
echo "Task id:$SLURM_ARRAY_TASK_ID"
echo "Variance experiments - Prepare settings for generating the populations in parallel"
echo "RNA: VRE & BE"

###First generate population
#Load required modules
module load mariadb
module load viennarna
module load java

#Start the Mysql Proxy
#ROOT_FOLDER=$PWD
#PROXYSCRIPT="$ROOT_FOLDER/proxy.sh"
#source $PROXYSCRIPT

#Prepare populations
./prepare_sample.sh $EXPERIMENT_ID