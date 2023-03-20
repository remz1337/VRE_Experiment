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
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=24G
##SBATCH --mem=160G
#SBATCH --output=slurm/analyze/slurm-%j.out
##SBATCH --error=slurm/analyze/slurm-%j.err

####### This process can take up to 200G of RAM when processing populations > 1M

EXPERIMENT_ID=$1
echo "Experiment id:$EXPERIMENT_ID"

echo "Slurm job id:$SLURM_JOB_ID"
#echo "Task id:$SLURM_ARRAY_TASK_ID"
echo "Variance experiments - Analyze populations with VRE program"
echo "RNA: VRE & BE"


#Make sure that the populations are generated beforehand (run generate script)

###Then run the measures
#Load required modules
#module load nixpkgs/16.09 gcc/6.4.0
module load gcc/11.3.0
module load mariadb

#Start the Mysql Proxy
#ROOT_FOLDER=$PWD
#PROXYSCRIPT="$ROOT_FOLDER/proxy.sh"
#source $PROXYSCRIPT

#To allow VRE to find libgfortran
#setrpaths.sh --path path/to/repo/VRE_Experiment/vre
#ROOT_FOLDER=$PWD
#setrpaths.sh --path $ROOT_FOLDER/vre

#Execute program
#./analyze.sh 1
#./analyze.sh 2
#./analyze.sh 3
#./analyze.sh 4
#./analyze.sh 5
#./analyze.sh 6
#./analyze.sh 7
#./analyze.sh 8
./analyze.sh $EXPERIMENT_ID
