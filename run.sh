#!/bin/bash
#PBS -q long
#PBS -l walltime=900:00:00
#PBS -l nodes=1:ppn=1
#PBS -k oe
#PBS -N Snakemake

module load python-3.8.2
cd /nlustre/users/graeme/pipeline-2020/
snakemake --cluster-config config/cluster.json --profile config/PBS-Torque-Profile 