#!/bin/sh
#SBATCH --job-name=drop	 	# Job name 			
#SBATCH --mail-type=ALL             	 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jpfeiffer@ufl.edu 	 # Where to send mail	
#SBATCH --nodes=1 	
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem-per-cpu=2g                 # mem per CPU
#SBATCH --time=60:00:00              # Time limit hrs:min:sec
#SBATCH --output=IBA_hog_%j.out     # Standard output and error log
#SBATCH --error=IBA_hog_%j.err		# Standard error log
#SBATCH --qos=page-b


module load python3/3.6.5

for M in T_D60_E15_A_FcC_* ; do python3 24c_flank_dropper.py $M FD22_$M B_platifrons 2 2;done