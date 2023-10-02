#!/bin/sh

#Move to exp/sample/walk directory
LOCALDIR=$(pwd)

echo "Local dir is $LOCALDIR"

cp vre.csv vre.csv.BAK

#First copy files to parent folder
for d in $(find $LOCALDIR/* -maxdepth 1 -type d)
do
  #Do something, the directory is accessible with $d:
  echo $d/vre.csv
  #sed '2d' $d/vre.csv >> $LOCALDIR/vre2.csv
  tail -n +3 $d/vre.csv >> $LOCALDIR/vre.csv
done

echo "Done!"