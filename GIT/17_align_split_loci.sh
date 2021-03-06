#!/bin/bash
#SBATCH --job-name=aaddlong
#SBATCH --output=mafft_array_%j.out
#SBATCH --error=maft_array_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jpfeiffer@ufl.edu
#SBATCH --time=10:00:00
#SBATCH --nodes=1-1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2G
#SBATCH --qos=page-b


module load mafft/7.245
	for X in L*.fa; do mafft --thread 8 $X > A_$X; done


