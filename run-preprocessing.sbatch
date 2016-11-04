#!/bin/bash
#SBATCH --job-name=pre-cc-1
#SBATCH --time=2-0:0:0
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --partition=dpetrov,hns
#################

export NXF_VER=0.22.4-SNAPSHOT
nextflow run rwtaylor/preprocessing -resume -c nextflow.config --samples samples.tsv -with-timeline -with-trace -with-dag preprocessing.png
