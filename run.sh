#!/bin/bash
#PBS -q long
#PBS -l walltime=900:00:00
#PBS -l nodes=1:ppn=14
#PBS -k oe
#PBS -o /nlustre/users/graeme/pipeline-2020/Logs/QSUB.log
#PBS -e /nlustre/users/graeme/pipeline-2020/Logs/QSUB.err
#PBS -N Snakemake Pipeline

module load python-3.8.2
cd /nlustre/users/graeme/pipeline-2020/
snakemake --cluster-config cluster.json --profile profile --report