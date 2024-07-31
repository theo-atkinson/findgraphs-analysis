#!/bin/bash
#PBS -P findgraphs-focused-scan
#PBS -q parallel20 
#PBS -l select=1:ncpus=20:mpiprocs=20:mem=90GB 
#PBS -j oe
#PBS -N focused-scan

cd [PATH TO FOCUSED SCAN R SCRIPT DIRECTORY]
cd $PBS_O_WORKDIR

# Load R v-4.2
source /app1/ebenv
module load R/4.2.2-foss-2022b

Rscript focused-scan.r
