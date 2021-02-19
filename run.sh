#!/bin/bash
#PBS -q long
#PBS -l walltime=900:00:00
#PBS -l nodes=1:ppn=1
#PBS -k oe
#PBS -o /nlustre/users/graeme/pipeline-msc/Logs/pipeline.log
#PBS -e /nlustre/users/graeme/pipeline-msc/Logs/pipeline.err
#PBS -N Snakemake Pipeline

module load python-3.8.2
cd /nlustre/users/graeme/pipeline-msc/
snakemake --cluster-config cluster.json --profile profile 
