#!/bin/sh

mkdir -p "slurm/vrna"

PHENOTYPE_LENGTH=30
PHENOTYPE_QUANTITY=100000

sbatch sbatch_vrna_generate_phenotypes.sh $PHENOTYPE_LENGTH $PHENOTYPE_QUANTITY
