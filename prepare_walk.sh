#!/bin/sh

#Define some global constants
ROOT_FOLDER=$PWD
EXPERIMENT_ID=$1
SAMPLE_ID=$2

good_walk_random_id="false"
while [ "$good_walk_random_id" == "false" ]; do
	#Generate random id
	WALK_RANDOM_ID=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1)

	WALK_RANDOM_ID_LOG="$SLURM_TMPDIR/$WALK_RANDOM_ID.log"

	$MYSQL_EXEC "insert into walk (fk_sample_id, random_id) values ($SAMPLE_ID, \"$WALK_RANDOM_ID\");" 2> "$WALK_RANDOM_ID_LOG"

	if [[ -z $(grep '[^[:space:]]' $WALK_RANDOM_ID_LOG) ]] ; then
		good_walk_random_id="true"
	fi
	#Cleanup
	rm $WALK_RANDOM_ID_LOG > /dev/null 2>&1
done

WALK_ID=$($MYSQL_EXEC "select ID from walk where fk_sample_id=$SAMPLE_ID and RANDOM_ID=\"$WALK_RANDOM_ID\";")

#abort if no walk id
if [ ${#WALK_ID} -eq 0 ]; then
	echo "Unable to get walk id for random id $WALK_RANDOM_ID."
	WALK_ID=$($MYSQL_EXEC "select ID from walk where fk_sample_id=$SAMPLE_ID and RANDOM_ID=\"$WALK_RANDOM_ID\";")
	echo "Try again. WALK ID:$WALK_ID"
fi

echo "Current walk id:$WALK_ID"

#Create subfolder FOR the walk
WALK_FOLDER="$SAMPLE_FOLDER/walk_$WALK_ID"
mkdir -p "$WALK_FOLDER"

TMP_WALK_FOLDER="$TMP_SAMPLE_FOLDER/walk_$WALK_ID"
mkdir -p "$TMP_WALK_FOLDER"

#Need to move into folder because ECJ outputs some files in execution directory
cd $TMP_WALK_FOLDER

VRE_POPULATION_FILE="$WALK_FOLDER/$VRE_POPULATION_FILENAME"
INITIAL_ECJ_POPULATION_FILE="$WALK_FOLDER/$INITIAL_ECJ_POPULATION_FILENAME"
TMP_INITIAL_ECJ_POPULATION_FILE="$TMP_WALK_FOLDER/$INITIAL_ECJ_POPULATION_FILENAME"
VALID_PHENOTYPES_FILE="$TMP_WALK_FOLDER/$VALID_PHENOTYPES_FILENAME"

#Update walk filename and sample size
$MYSQL_EXEC "update walk set initial_population_size = $POPULATION_SIZE where ID = $WALK_ID;"

echo "Updated walk id initial population size to:$POPULATION_SIZE"

#Setting up ahead so it's available for PREPARE_BNK
ECJ_RUN_ARGUMENTS="$ECJ_BASE_ARGUMENTS"

if [ "$EXPERIMENT_TYPE" == "vre" ]; then

	#Make sure to keep the offspring
	BASELINE_KEEP_BEST=false
	ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p vre.baselinemutation.keepbest=$BASELINE_KEEP_BEST"

	#Setup inital population (generation 0)
	if [ "$ALGORITHM_NAME" == "$INVERSERNA_NAME" ]; then
		echo "Inverse RNA"
		population_counter=0
		PHENOTYPE_QTY=1

		RANDOM_PHENOTYPE=$($MYSQL_EXEC "BEGIN; SET @selected_phenotype := '0'; UPDATE available_phenotypes SET used = 1, phenotype = (SELECT @selected_phenotype := phenotype) WHERE length = $GENOTYPE_LENGTH and used =0 LIMIT 1; SELECT @selected_phenotype; COMMIT;")
		echo "Selected phenotype:$RANDOM_PHENOTYPE"
		source $VRNA_PHENOTYPE_REGENERATOR $RANDOM_PHENOTYPE $GENOTYPE_LENGTH $INVERSE_REPEAT

		#just in case the inverse mapping fails to find at least one valid genotype
		if [ $population_counter -le 0 ]; then
			echo "For some reason the population counter is still 0... Generating a new phenotype"
			source $VRNA_PHENOTYPE_GENERATOR $PHENOTYPE_QTY $GENOTYPE_LENGTH $INVERSE_REPEAT
		fi

		ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p eval.problem.target-phenotype=$RANDOM_PHENOTYPE"

		#Prepare population file in ECJ format
		$VRE_ECJ_BUILDER -i $VRE_POPULATION_FILE -c > $TMP_INITIAL_ECJ_POPULATION_FILE

#### Not used anymore
#		#Number of lines to remove from each ECJ run CSV output
#		HEADER_SIZE=$population_counter
#		#((HEADER_SIZE++))
#		#We need to do +2 here because when using the tail function to skip only a few lines at the beginning, the counter starts at -1
#		#So +1 to negate that, then another +1 for the actual csv header in the file
#		((HEADER_SIZE+=2))

		#Update sample size
		$MYSQL_EXEC "update walk set initial_population_size = $population_counter where ID = $WALK_ID;"
#### I don't think I use this anymore... Should clean up the hardcoded names
#	elif [ "$ALGORITHM_NAME" == "$BNK_OSC_VRE_NAME" ]; then
#		echo "BNK Oscillator VRE"
#		#Create initial pop from an available genotype
#		source $PREPARE_BNK
	else
		echo "VRE"

		#Create starting population (generation 0) from ECJ
		ECJ_SETUP_ARGUMENTS="$ECJ_BASE_ARGUMENTS -p stat.VRE-CSV=$VRE_POPULATION_FILE -p generations=1 -p pop.subpop.0.size=$POPULATION_SIZE -p finish=ec.vre.VREFinisher"
		java $ECJ_SETUP_ARGUMENTS > /dev/null 2>&1

#####Not used anymore
#		#Number of lines to remove from each ECJ run CSV output
#		HEADER_SIZE=$POPULATION_SIZE
#		#((HEADER_SIZE++))
#		((HEADER_SIZE+=2))
	fi

	#Copy initial pop file to network storage
	cp $TMP_INITIAL_ECJ_POPULATION_FILE $INITIAL_ECJ_POPULATION_FILE
	ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p pop.file=$INITIAL_ECJ_POPULATION_FILE -p breed=ec.vre.VREBreeder -p pop.subpop.0.species.pipe=ec.vector.breed.VectorNeighborhoodMutationPipeline"
	ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p quit-on-run-complete=false"

	neighborhood_it=1
	while [ $neighborhood_it -le $NEIGHBORHOOD ]; do
		#Create subfolder for the neighborhood
		NEIGHBORHOOD_FOLDER="$WALK_FOLDER/neighborhood_$neighborhood_it"
		mkdir -p "$NEIGHBORHOOD_FOLDER"

		TMP_NEIGHBORHOOD_FOLDER="$TMP_WALK_FOLDER/neighborhood_$neighborhood_it"
		mkdir -p "$TMP_NEIGHBORHOOD_FOLDER"
		cd $TMP_NEIGHBORHOOD_FOLDER

		NEIGHBORHOOD_OUTPUT_FILE="$NEIGHBORHOOD_FOLDER/$VRE_POPULATION_FILENAME"
		MATH_HELPER_FILE="$TMP_NEIGHBORHOOD_FOLDER/$MATH_HELPER_FILENAME"

		#Do not waste time oversampling small neighborhoods
		rm $MATH_HELPER_FILE > /dev/null 2>&1
		MATH_HELP_ARGS="-n $neighborhood_it -l $GENOTYPE_LENGTH -p $PROB_NAME"
		if [ "$PROB_NAME" == "RNA" ]; then
			MATH_HELP_ARGS="$MATH_HELP_ARGS -a $RNA_ALPHABET_SIZE"
		elif [ "$PROB_NAME" == "BNK" ]; then
			MATH_HELP_ARGS="$MATH_HELP_ARGS -g $BNK_GATES -i $BNK_INPUTS"
		fi

		$VRE_MATH_HELPER $MATH_HELP_ARGS > $MATH_HELPER_FILE
		MATH_HELPER_OUTPUT=$(head -n 1 $MATH_HELPER_FILE)
		SEARCH_SPACE=$(tail -n 1 $MATH_HELPER_FILE)
		
		#echo "SEARCH_SPACE:$SEARCH_SPACE"

		MAX_NEIGHBORHOOD_GENERATIONS=$(echo $MATH_HELPER_OUTPUT | cut -f1 -d\|)
		#SAMPLING_RATIO=$(echo $MATH_HELPER_OUTPUT | cut -f2 -d\|)

		RUN_NEIGHBORHOOD_GENERATIONS=$GENERATIONS
		### Simply do either the entire space if it's less than the predefined max generations
#		if [ "$GENERATIONS" -gt "$MAX_NEIGHBORHOOD_GENERATIONS" ]; then
#			RUN_NEIGHBORHOOD_GENERATIONS=$MAX_NEIGHBORHOOD_GENERATIONS
		if [ "$GENERATIONS" -gt "$SEARCH_SPACE" ]; then
			RUN_NEIGHBORHOOD_GENERATIONS=$SEARCH_SPACE
		fi

		#Explore each neighborhood, from previously generated initial population
		ECJ_RUN_ARGUMENTS_FINAL="$ECJ_RUN_ARGUMENTS -p generations=$RUN_NEIGHBORHOOD_GENERATIONS -p stat.VRE-CSV=$NEIGHBORHOOD_OUTPUT_FILE -p vre.baselinemutation.neighborhood=$neighborhood_it"

		#Queue job
		$MYSQL_EXEC "replace into vre_neighborhood_job (distance, folder, population_file, java_args, generated_neighbors, fk_walk_id) values ($neighborhood_it, \"$NEIGHBORHOOD_FOLDER\", \"$NEIGHBORHOOD_OUTPUT_FILE\", \"$ECJ_RUN_ARGUMENTS_FINAL\", $RUN_NEIGHBORHOOD_GENERATIONS, $WALK_ID);"

		((neighborhood_it++))
	done

else
	ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p generations=$GENERATIONS -p pop.subpop.0.size=$POPULATION_SIZE -p stat.VRE-CSV=$VRE_POPULATION_FILE"

	#Leverage same starting points (RNA)
	
	if [ ! -z "$BNK_SOURCE_EXP_ID" ] && [ $BNK_SOURCE_EXP_ID -gt 0 ]; then
		if [ "$PROB_NAME" == "RNA" ]; then
			#if a source experiment is provided, then use the same starting points
			source $PREPARE_RNA

#			ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p pop.file=$INITIAL_ECJ_POPULATION_FILE"

#			SOURCE_WALK_FILE="$WALK_FOLDER/source_walk.txt"
#			#Log the source walk ID from the BE experiment used to initialize the BR walk
#			echo "$SOURCE_WALK_ID" > $SOURCE_WALK_FILE
		elif [ "$PROB_NAME" == "BNK" ]; then
			source $PREPARE_BNK

#			#copy to persistent storage
#			cp "$TMP_INITIAL_ECJ_POPULATION_FILE" "$INITIAL_ECJ_POPULATION_FILE"
#			ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p pop.file=$INITIAL_ECJ_POPULATION_FILE"
			ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p pop.subpop.0.species.ind.n=$BNK_GATES -p pop.subpop.0.species.ind.k=$BNK_INPUTS"

#			SOURCE_WALK_FILE="$WALK_FOLDER/source_walk.txt"
#			#Log the source walk ID from the BE experiment used to initialize the BR walk
#			echo "$SOURCE_WALK_ID" > $SOURCE_WALK_FILE	
		else
			source $PREPARE_BNK

#			#copy to persistent storage
#			cp "$TMP_INITIAL_ECJ_POPULATION_FILE" "$INITIAL_ECJ_POPULATION_FILE"
#			ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p pop.file=$INITIAL_ECJ_POPULATION_FILE"

#			SOURCE_WALK_FILE="$WALK_FOLDER/source_walk.txt"
#			#Log the source walk ID from the BE experiment used to initialize the BR walk
#			echo "$SOURCE_WALK_ID" > $SOURCE_WALK_FILE
		fi
		
		ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p pop.file=$INITIAL_ECJ_POPULATION_FILE"
		
		SOURCE_WALK_FILE="$WALK_FOLDER/source_walk.txt"
		#Log the source walk ID from the BE experiment used to initialize the BR walk
		echo "$SOURCE_WALK_ID" > $SOURCE_WALK_FILE
	fi

	#Run the experiment
	#If baseline evolvability, need to generate the target phenotype randomly
	if [ "$EXPERIMENT_TYPE" == "be" ]; then
		echo "BE"
		ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p quit-on-run-complete=true"

		BASELINE_KEEP_BEST=true
		PHENOTYPE_QTY=1
		if [ "$PROB_NAME" == "RNA" ]; then
			#First, try to pull one from the preprocessed VRNA DB
			TARGET_PHENOTYPE=$($MYSQL_EXEC "BEGIN; SET @selected_phenotype := '0'; UPDATE available_phenotypes SET used = 1, phenotype = (SELECT @selected_phenotype := phenotype) WHERE length = $GENOTYPE_LENGTH and used =0 LIMIT 1; SELECT @selected_phenotype; COMMIT;")
			echo "Target phenotype from DB:$TARGET_PHENOTYPE"
			if [ ${#TARGET_PHENOTYPE} -ne $GENOTYPE_LENGTH ]; then
			#generate one
				echo "For some reason the target phenotype doesn't fit the genotype length... Generating a new phenotype"
				source $VRNA_PHENOTYPE_GENERATOR $PHENOTYPE_QTY $GENOTYPE_LENGTH $INVERSE_REPEAT > $VALID_PHENOTYPES_FILE
				TARGET_PHENOTYPE=$(head -n 1 $VALID_PHENOTYPES_FILE)
			fi

		elif [ "$PROB_NAME" == "BNK" ]; then
			echo "BNK Oscillator/Walk"
			ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p pop.subpop.0.species.ind.n=$BNK_GATES -p pop.subpop.0.species.ind.k=$BNK_INPUTS"
		fi

		ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p eval.problem.target-phenotype=$TARGET_PHENOTYPE"
	else
		echo "BR"

		ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p quit-on-run-complete=false"
		
		BASELINE_KEEP_BEST=false

#		if [ "$PROB_NAME" == "BNK" ]; then
#			source $PREPARE_BNK

#			#copy to persistent storage
#			cp "$TMP_INITIAL_ECJ_POPULATION_FILE" "$INITIAL_ECJ_POPULATION_FILE"
#			ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p pop.file=$INITIAL_ECJ_POPULATION_FILE -p pop.subpop.0.species.ind.n=$BNK_GATES -p pop.subpop.0.species.ind.k=$BNK_INPUTS"

#			SOURCE_WALK_FILE="$WALK_FOLDER/source_walk.txt"
#			#Log the source walk ID from the BE experiment used to initialize the BR walk
#			echo "$SOURCE_WALK_ID" > $SOURCE_WALK_FILE
#		fi
	fi

	ECJ_RUN_ARGUMENTS="$ECJ_RUN_ARGUMENTS -p vre.baselinemutation.keepbest=$BASELINE_KEEP_BEST -p breed=ec.simple.SimpleBreeder"

	#Queue job
	$MYSQL_EXEC "replace into baseline_job (folder, population_file, java_args, target_phenotype, fk_walk_id) values (\"$WALK_FOLDER\", \"$VRE_POPULATION_FILE\", \"$ECJ_RUN_ARGUMENTS\", \"$TARGET_PHENOTYPE\", $WALK_ID);"
fi
