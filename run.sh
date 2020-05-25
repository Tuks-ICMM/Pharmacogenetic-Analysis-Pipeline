#!/bin/bash
#PBS -q normal
#PBS -l walltime=05:00:00
#PBS -l nodes=1:ppn=24
#PBS -k oe
#PBS -o /nlustre/users/graeme/pipeline-2020/Logs/QSUB.log
#PBS -e /nlustre/users/graeme/pipeline-2020/Logs/QSUB.err
#PBS -N Snakemake Pipeline

module load python-3.6.6
snakemake --cluster-config cluster.json