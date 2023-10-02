#!/bin/sh
#SBATCH --job-name=Plot
#SBATCH --time=01-08:00:00
#SBATCH --account=slurm-account
#SBATCH --mail-user=myemail@email.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=4
##SBATCH --mem-per-cpu=12G
#SBATCH --mem=160G
##SBATCH --array=2-4
#SBATCH --output=slurm/plotter/slurm-%j.out
##SBATCH --error=slurm/plotter/slurm-%j.err

### For HUGE experiments (such as the BNK(4,4) multi-modality with 10M individuals), run for 5 days, with 120G RAM

EXPERIMENT_ID=$1
echo "Experiment id:$EXPERIMENT_ID"

echo "Slurm job id:$SLURM_JOB_ID"
#echo "Task id:$SLURM_ARRAY_TASK_ID"
echo "VRE Plotter"
echo "Plot graphs"

#Load required modules
module load scipy-stack

#Start the Mysql Proxy
#ROOT_FOLDER=$PWD
#PROXYSCRIPT="$ROOT_FOLDER/proxy.sh"
#source $PROXYSCRIPT

#Use imageio to create gifs
pip3 install imageio
pip3 install mysql-connector-python
pip3 install graphviz
#pip3 install transitions
#pip3 install sqlalchemy
pip3 install 'sqlalchemy<2.0'
pip3 install scipy

#Execute program
./plot.sh $EXPERIMENT_ID
