#!/bin/bash


#SBATCH --job-name=QTRIM
#SBATCH -o %A_%a.IBA_P3.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=pfeifferj@si.edu
#SBATCH --mem-per-cpu=4000mb
#SBATCH -t 01:00:00
#SBATCH -n 1
#SBATCH -c 1

#SBATCH --qos=page-b
#SBATCH --array=1-200

module load cutadapt trim_galore

INFILE1=$(head -n $SLURM_ARRAY_TASK_ID RAWLIST | tail -n1|cut -f1)
INFILE2=$(head -n $SLURM_ARRAY_TASK_ID RAWLIST | tail -n1|cut -f2)


trim_galore -q 20 --phred33 --length 30 --dont_gzip --paired $INFILE1 $INFILE2
