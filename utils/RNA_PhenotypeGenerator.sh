#!/bin/sh

#These variables need to be SET by the calling script:
#ROOT_FOLDER, SAMPLE_FOLDER, VRE_CSV_FILE
#Assume we should already be in the right subfolder (sample/walk) to create folders for each phenotype
#START_DIR=$PWD
#UTILS_DIR=$(dirname $(readlink -f $0))
UTILS_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
ROOT_DIR=$(dirname "$UTILS_DIR")


#FOR VRNA
QUANTITY=$1
LENGTH=$2
INVERSE_REPEAT=$3

#Define some constant filenames
STRUCTURE_FILENAME="structures.txt"
INVERSE_FILENAME="inverse.out"
REFOLD_FILENAME="refold.fa"
REFOLDED_FILENAME="rna_sequence.out"
#SET VRE_CSV_FILENAME=rnavre_util.csv
#SET VALID_PHENOTYPES_FILENAME=valid_phenotypes.txt

VRE_PHENOTYPE_GENERATOR="$UTILS_DIR/VRE_PhenotypeGenerator"
VRNA_PARSER="$ROOT_DIR/vrna/VRNA_Parser"
#SET VRE_CSV_FILE=!SAMPLE_FOLDER!\%VRE_CSV_FILENAME%
#SET VALID_PHENOTYPES_FILE=!SAMPLE_FOLDER!\%VALID_PHENOTYPES_FILENAME%

echo Genotype,Fitness,Phenotype > $VRE_POPULATION_FILE

structure_counter=0
population_counter=0
last_population_size=0
while [ $population_counter -lt $QUANTITY ]; do
	#create subdir for the structure
	#STRUCT_FOLDER="$WALK_FOLDER/structures/struct_$structure_counter"
	STRUCT_FOLDER="$TMP_WALK_FOLDER/structures/struct_$structure_counter"
	mkdir -p "$STRUCT_FOLDER"

	#Define local files
	STRUCTURE_FILE="$STRUCT_FOLDER/$STRUCTURE_FILENAME"
	INVERSE_FILE="$STRUCT_FOLDER/$INVERSE_FILENAME"
	REFOLD_FILE="$STRUCT_FOLDER/$REFOLD_FILENAME"
	REFOLDED_FILE="$STRUCT_FOLDER/$REFOLDED_FILENAME"
	#VRE_CSV_FILE=!STRUCT_FOLDER!\!VRE_CSV_FILENAME!

	#call c++ exe to generate phenotype (secondary structures)
	$VRE_PHENOTYPE_GENERATOR -p RNA -l $LENGTH -q 1 > $STRUCTURE_FILE

	#read secondary structure from file
	structure=$(head -n 1 $STRUCTURE_FILE)
	#echo "struct $structure"

	#execute rnainverse
	line_count=0

	while [ $line_count -lt $INVERSE_REPEAT ]; do
		#RNAInverse -Fmp -f 0.99 -R%INVERSE_REPEAT% < %STRUCTURE_FILE% > %INVERSE_FILE%
		RNAinverse -R$INVERSE_REPEAT < $STRUCTURE_FILE > $INVERSE_FILE

		line_count=$(< "$INVERSE_FILE" wc -l)
	done

	#First, Refold to validate the generate RNA sequence from the inverse folder
	$VRNA_PARSER -i $INVERSE_FILE -r --target="$structure" > $REFOLD_FILE

	##execute rnafold
	#RNAfold --noPS < $REFOLD_FILE > $REFOLDED_FILE
	RNAsubopt -s < $REFOLD_FILE > $REFOLDED_FILE

	$VRNA_PARSER -i $INVERSE_FILE -f $REFOLDED_FILE --target="$structure" >> $VRE_POPULATION_FILE

	#read number of lines in CSV and continue loop until population size >= quantity
	population_counter=$(< "$VRE_POPULATION_FILE" wc -l)
	#start pop at -1 since we need to exclude the first line in the csv file
	((population_counter--))

	#Need to review this... Does not mean it is a unique phenotype because it worked...
	if [ $population_counter -gt $last_population_size ]; then
		#echo $structure>>"$VALID_PHENOTYPES_FILE"
		echo $structure
	fi
	
	#cleanup
	rm -rf $STRUCT_FOLDER

	last_population_size=$population_counter
	((structure_counter++))
done
