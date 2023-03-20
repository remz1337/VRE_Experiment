#!/bin/sh

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")
VRNA_Evaluator="$SCRIPTPATH/VRNA_Parser"

#read rna from stdin
inputfile="$1"
target="$2"
outputfile="$3"

SEQUENCE_FILE="vrna_sequence.out"

#execute rnafold
#RNAfold --noPS -p -d2 --noDP --MEA  < "$inputfile" > "$SEQUENCE_FILE"
RNAsubopt -s < "$inputfile" > "$SEQUENCE_FILE"

#call c++ exe
$VRNA_Evaluator -e -f "$SEQUENCE_FILE" -t "$target" > "$outputfile"
