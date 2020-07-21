#!/bin/bash
#SBATCH --job-name=TRIMAL
#SBATCH --output=trimal_%j.out
#SBATCH --error=trimal_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jpfeiffer@ufl.edu
#SBATCH --time=12:00:00
#SBATCH --nodes=1-1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2G
#SBATCH --qos=page-b
#SBATCH --array 1-539

module load trimal/1.2

input_aln=`ls *.fas | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

out=`echo $input_aln`.trim_fas
htmlout=`echo $input_aln`.trim_html

trimal -in $input_aln -out $out -htmlout $htmlout -gt 0.5