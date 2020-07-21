#!/bin/bash
#SBATCH --job-name=gene_tree
#SBATCH --output=raxml_bash_%j.out
#SBATCH --error=raxml_bash_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jpfeiffer@ufl.edu
#SBATCH --time=24:00:00
#SBATCH --nodes=1-1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=2G
#SBATCH --qos=page-b
#SBATCH --array 1-579%10

module load iq-tree/1.6.12

input_aln=`ls *.fas | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
out=`echo $input_aln`.tre

iqtree -s $input_aln -bb 1000