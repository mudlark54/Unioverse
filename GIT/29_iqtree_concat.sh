#!/bin/bash
#SBATCH --job-name=cat
#SBATCH --output=iqtree_%j.out
#SBATCH --error=iqtree_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jpfeiffer@ufl.edu
#SBATCH --time=60:00:00
#SBATCH --nodes=1-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --qos=page-b

module load iq-tree/1.6.12

iqtree -nt AUTO -s FcC_supermatrix.phy -spp FcC_supermatrix_partition.txt -m MF+MERGE -rcluster 10 -bb 1000 -safe