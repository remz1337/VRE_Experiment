#!/bin/sh

#Define some constant filenames
STRUCTURE_FILENAME="structures.txt"
INVERSE_FILENAME="inverse.out"
REFOLD_FILENAME="refold.fa"
REFOLDED_FILENAME="rna_sequence.out"
VALID_PHENOTYPES_FILENAME="valid_phenotypes.txt"
JAVA_LOG_FILENAME="javaoutput.log"
RUN_LOG_FILENAME="output.log"
INITIAL_ECJ_POPULATION_FILENAME="population.in"
SQLITE_ERROR_FILENAME="sqlite_errors.err"
MATH_HELPER_FILENAME="vre_math_helper.txt"
SEQUENCE_FILE="vrna_sequence.out"
VRNA_PARSER_DISTANCE="distances.out"
VRNA_PARSER_INPUT="rna_sequence.fa"


#Clean up temporary files
rm $STRUCTURE_FILENAME > /dev/null 2>&1
rm $INVERSE_FILENAME > /dev/null 2>&1
rm $REFOLD_FILENAME > /dev/null 2>&1
rm $REFOLDED_FILENAME > /dev/null 2>&1
rm $VALID_PHENOTYPES_FILENAME > /dev/null 2>&1
rm $JAVA_LOG_FILENAME > /dev/null 2>&1
rm $RUN_LOG_FILENAME > /dev/null 2>&1
rm $INITIAL_ECJ_POPULATION_FILENAME > /dev/null 2>&1
rm $SQLITE_ERROR_FILENAME > /dev/null 2>&1
rm $MATH_HELPER_FILENAME > /dev/null 2>&1
rm $SEQUENCE_FILE > /dev/null 2>&1

rm $VRNA_PARSER_DISTANCE > /dev/null 2>&1
rm $VRNA_PARSER_INPUT > /dev/null 2>&1
