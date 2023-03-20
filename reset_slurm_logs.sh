#!/bin/sh

#Clean up slurm logs
rm -rf slurm/prepare
rm -rf slurm/prepare_sample
rm -rf slurm/generate
rm -rf slurm/generate_vre
rm -rf slurm/generate_baseline
rm -rf slurm/analyze
rm -rf slurm/plotter
rm -rf slurm/compress
#rm -rf slurm/proxy

#Create slurm dirs
mkdir -p "slurm/prepare"
mkdir -p "slurm/prepare_sample"
mkdir -p "slurm/generate"
mkdir -p "slurm/generate_vre"
mkdir -p "slurm/generate_baseline"
mkdir -p "slurm/analyze"
mkdir -p "slurm/plotter"
mkdir -p "slurm/compress"
mkdir -p "slurm/proxy"
