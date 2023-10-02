#!/bin/sh
##SBATCH --job-name=SINGLETON
#SBATCH --time=00-00:30:00
#SBATCH --account=slurm-account
#SBATCH --mail-user=myemail@email.com
##SBATCH --mail-type=BEGIN
##SBATCH --mail-type=END
##SBATCH --mail-type=FAIL
##SBATCH --mail-type=REQUEUE
##SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=slurm/prepare/slurm-%j.out
##SBATCH --error=slurm/prepare/slurm-%j.err

EXPERIMENT_ID=$1
echo "Experiment id:$EXPERIMENT_ID"

echo "Slurm job id:$SLURM_JOB_ID"
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
./prepare.sh $EXPERIMENT_ID