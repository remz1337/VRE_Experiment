#!/bin/sh
#SBATCH --job-name=Generate_VRNA_Phenotypes
#SBATCH --time=00-16:00:00
#SBATCH --account=slurm-account
#SBATCH --mail-user=myemail@email.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
##SBATCH --array=0-3
#SBATCH --output=slurm/vrna/slurm-%j.out
##SBATCH --error=slurm/vrna/slurm-%j.err

echo "Slurm job id:$SLURM_JOB_ID"
#echo "Task id:$SLURM_ARRAY_TASK_ID"
echo "Generate RNA phenotypes"

###First generate population
#Load required modules
module load viennarna
module load java

RNA_LENGTH=$1
QUANTITY=$2

#Generate RNA phenotypes
./vrna_generate_phenotypes.sh $RNA_LENGTH $QUANTITY $SLURM_JOB_ID
