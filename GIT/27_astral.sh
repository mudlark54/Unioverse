#!/bin/bash
#SBATCH --job-name=astral
#SBATCH --output=iqtree_%j.out
#SBATCH --error=iqtree_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jpfeiffer@ufl.edu
#SBATCH --time=40:00:00
#SBATCH --nodes=1-1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=3G
#SBATCH --qos=page-b


module load astral/5.6.1

astral -i 1_alltrees_sans10 -o 1_alltrees_sans10.out 2> 1_alltrees_sans10.astral.log

