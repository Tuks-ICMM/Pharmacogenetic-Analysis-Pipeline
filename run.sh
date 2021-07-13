#!/bin/bash
#PBS -q long
#PBS -l walltime=900:00:00
#PBS -l nodes=1:ppn=1
#PBS -k oe
#PBS -o /nlustre/users/graeme/pipeline-2020/logs/pipeline.log
#PBS -e /nlustre/users/graeme/pipeline-2020/logs/pipeline.err
#PBS -N Snakemake Pipeline

module load anaconda3-2020.07
cd /nlustre/users/graeme/pipeline-2020/
snakemake --cluster-config cluster.json --profile profile 