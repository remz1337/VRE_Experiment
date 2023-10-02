#!/bin/sh
##SBATCH --job-name=SINGLETON
#SBATCH --time=00-01:00:00
#SBATCH --account=slurm-account
#SBATCH --mail-user=myemail@email.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=slurm/generate/slurm-%j.out
##SBATCH --error=slurm/generate/slurm-%j.err

EXPERIMENT_ID=$1
echo "Experiment id:$EXPERIMENT_ID"

echo "Slurm job id:$SLURM_JOB_ID"
echo "Variance experiments - Generate populations"
echo "RNA: VRE & BE"

###First generate population
#Load required modules
module load viennarna
module load java
module load mariadb

#Start the Mysql Proxy
ROOT_FOLDER=$PWD
#PROXYSCRIPT="$ROOT_FOLDER/proxy.sh"
#source $PROXYSCRIPT

LOAD_VARIABLES="$ROOT_FOLDER/load_variables.sh"
source $LOAD_VARIABLES

#ROOT_FOLDER=$PWD
#GENERATE_SCRIPT="$ROOT_FOLDER/generate.sh"

source $PREPARE_LOAD $EXPERIMENT_ID
#Dynamically adjust memory needs, based on population size and number of generation
MEM_PER_CPU="2G"
JOB_TIME="00-06:00:00"
CPUS=6
EXP_SCALE=$((POPULATION_SIZE * GENERATIONS))
if [ $EXP_SCALE -ge 1000000 ]; then 
	MEM_PER_CPU="24G"
	JOB_TIME="01-00:00:00"
	CPUS=6
elif [ $EXP_SCALE -ge 100000 ]; then 
	MEM_PER_CPU="12G"
	JOB_TIME="00-12:00:00"
	CPUS=8
elif [ $EXP_SCALE -ge 10000 ]; then 
	MEM_PER_CPU="6G"
	JOB_TIME="00-08:00:00"
	CPUS=12
else
	MEM_PER_CPU="4G"
	JOB_TIME="00-04:00:00"
	CPUS=16
fi

source $GENERATE_SCRIPT $EXPERIMENT_ID

#Wait for all java job arrays to complete, then run the analysis job
#ANALYZE_JOB_ID=$(sbatch --parsable --dependency=afterok:$GENERATE_JOB_ID sbatch_analyze.sh $EXPERIMENT_ID)
ANALYZE_JOB_ID=$(sbatch --mem-per-cpu=$MEM_PER_CPU --time=$JOB_TIME --cpus-per-task=$CPUS --parsable --dependency=afterany:$GENERATE_JOB_ID sbatch_analyze.sh $EXPERIMENT_ID)

#Moved into analyze script
#Then create the plots
#sbatch --dependency=afterok:$ANALYZE_JOB_ID sbatch_plotter.sh $EXPERIMENT_ID
#PLOTTER_JOB_ID=$(sbatch --parsable --dependency=afterok:$ANALYZE_JOB_ID sbatch_plotter.sh $EXPERIMENT_ID)

#Finally, compress and delete files
#COMPRESS_JOB_ID=$(sbatch --parsable --dependency=afterok:$PLOTTER_JOB_ID sbatch_compress.sh $EXPERIMENT_ID)
