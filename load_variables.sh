#!/bin/sh

#---------------- Variables

#Define some global constants
ROOT_FOLDER=$PWD
INVERSERNA_NAME="IMAP"
####Shouldn't use these BNK hardcoded names anymore... clean up
#BNK_OSC_BE_NAME="ECJ_BNK_BE"
#BNK_WALK_BE_DEMO_NAME="ECJ_BNK_BE_DEMO"
#BNK_OSC_VRE_NAME="ECJ_BNK_VRE"
SYSTEM_INTEGER_LIMIT=2147483647
#MAX_GENERATIONS_WITH_BUFFER=2045222520
#Following method from paper "Determining sample size...", max shouldn't be above ~16k. 20k should be enough to cover for duplicate entries
MAX_GENERATIONS_WITH_BUFFER=20000

#Common filenames
#INITIAL_VRE_POPULATION_FILENAME="rnavre.csv"
#FINAL_POPULATION_FILENAME="vre_pop.csv"
VRE_POPULATION_FILENAME="vre.csv"
VALID_PHENOTYPES_FILENAME="valid_phenotypes.txt"
INITIAL_ECJ_POPULATION_FILENAME="population.in"
SQLITE_ERROR_FILENAME="sqlite_errors.err"
MATH_HELPER_FILENAME="vre_math_helper.txt"
#MYSQL_CREDENTIALS_FILENAME="mysql_creds.env"
#The config file for MariaDB (MySQL) should be like this:
#[client] 
#host=<your-db-host>
#port=<your-db-port>
#socket=<your-db-socket>
#database=<your-db-name>
#user=<your-db-user>
#password=<your-db-password>


#Tools
VRE="$ROOT_FOLDER/vre/VRE"
VRE_ECJ_BUILDER="$ROOT_FOLDER/vre/VRE_ECJ_Builder"
VRNA_PHENOTYPE_GENERATOR="$ROOT_FOLDER/utils/RNA_PhenotypeGenerator.sh"
VRNA_PHENOTYPE_REGENERATOR="$ROOT_FOLDER/utils/RNA_PhenotypeRegenerator.sh"
VRE_MATH_HELPER="$ROOT_FOLDER/utils/VRE_MathHelper"
VRE_PHENOTYPE_GENERATOR="$ROOT_FOLDER/utils/VRE_PhenotypeGenerator"
PREPARE_SAMPLE="$ROOT_FOLDER/prepare_sample.sh"
PREPARE_WALK="$ROOT_FOLDER/prepare_walk.sh"
PREPARE_BNK="$ROOT_FOLDER/prepare_bnk.sh"
PREPARE_RNA="$ROOT_FOLDER/prepare_rna.sh"
GENERATE_SCRIPT="$ROOT_FOLDER/generate.sh"
MYSQL_EXEC="$ROOT_FOLDER/mysql_exec.sh"
MYSQL_EXEC_JOB="$ROOT_FOLDER/mysql_exec_job.sh"

PREPARE_LOAD="$ROOT_FOLDER/prepare_load.sh"
SBATCH_PREPARE_SAMPLE="$ROOT_FOLDER/sbatch_prepare_sample.sh"
SBATCH_GENERATE="$ROOT_FOLDER/sbatch_generate.sh"

#Common filenames
MYSQL_CREDENTIALS_FILENAME="mysql_creds.env"

#Executables
#VRE="$ROOT_FOLDER/vre/VRE"
#MYSQL_EXEC="$ROOT_FOLDER/mysql_exec.sh"
#PREPARE_LOAD="$ROOT_FOLDER/prepare_load.sh"

#Load MySQL credentials
MYSQL_CREDENTIALS_FILE="$ROOT_FOLDER/$MYSQL_CREDENTIALS_FILENAME"




#ECJ Variables
RNG_SEED=time
BASELINE_KEEP_BEST=false

#RNAInverse
INVERSE_REPEAT=50