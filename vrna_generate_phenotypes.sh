#!/bin/sh

#---------------- Variables

#Define some global constants
ROOT_FOLDER=$PWD
JOB_ID=$3

#Common filenames
#INITIAL_VRE_POPULATION_FILENAME="rnavre.csv"
VRE_POPULATION_FILENAME="vre.csv"
JAVA_LOG_FILENAME="javaoutput.log"
RUN_LOG_FILENAME="output.log"
SQLITE_ERROR_FILENAME="sqlite_errors.err"
EXTRACTED_PHENOTYPES_FILENAME="extracted_phenotypes.csv"

#Tools
VRNA_DATABASE="$ROOT_FOLDER/vrna/vrna.db3"
VRNA_EXTRACTOR="$ROOT_FOLDER/vrna/VRNA_Extractor"

#RNA Length
RNA_LENGTH=$1

#Number of phenotypes to generate (make this high, job will likely timeout when resources are exhausted)
PHENOTYPE_QTY=$2

VRNA_FOLDER="$ROOT_FOLDER/vrna/job_$JOB_ID"
mkdir -p "$VRNA_FOLDER"
cd $VRNA_FOLDER

VRE_POPULATION_FILE="$VRNA_FOLDER/$VRE_POPULATION_FILENAME"
JAVA_LOG_FILE="$VRNA_FOLDER/$JAVA_LOG_FILENAME"
EXTRACTED_PHENOTYPES_FILE="$VRNA_FOLDER/$EXTRACTED_PHENOTYPES_FILENAME"

ECJ_ARGS="-jar $ROOT_FOLDER/jars/ecj_rna.jar -file $ROOT_FOLDER/params/rna/rna.simple.params -p seed.0=time -p pop.subpop.0.species.genome-size=$RNA_LENGTH -p pop.subpop.0.size=$PHENOTYPE_QTY -p generations=1 -p stat.VRE-CSV=$VRE_POPULATION_FILE -p breed=ec.simple.SimpleBreeder -p pop.subpop.0.species.pipe=ec.vector.breed.VectorNeighborhoodMutationPipeline"

#Build initial random genotypes and fold them to get the phenotypes
java $ECJ_ARGS &> "$JAVA_LOG_FILE"

#Extract only the phenotypes (no duplication) for easy insert in the DB
$VRNA_EXTRACTOR -i "$VRE_POPULATION_FILE" > "$EXTRACTED_PHENOTYPES_FILE"

#sqlite3 "$VRNA_DATABASE" ".mode csv; .import $EXTRACTED_PHENOTYPES_FILE available_phenotypes"
#sqlite3 -separator ',' "$VRNA_DATABASE" ".import $EXTRACTED_PHENOTYPES_FILE available_phenotypes"
#Import the CSV files manually in MySQL


#Insert in DB
#while IFS= read -r line
#do
#	if [ ${#line} -eq $RNA_LENGTH ]; then 
#		sqlite3 "$VRNA_DATABASE" "insert or ignore into available_phenotypes (phenotype, length) values (\"$line\", $RNA_LENGTH);"
#	fi
#done < "$EXTRACTED_PHENOTYPES_FILE"


cd $ROOT_FOLDER
echo Done!