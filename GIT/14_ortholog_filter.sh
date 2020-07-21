#!/bin/sh
#SBATCH --job-name=otholog_filter	 	# Job name 			
#SBATCH --mail-type=ALL             	 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jpfeiffer@ufl.edu 	 # Where to send mail	
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=2                   # Run a single task	
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem-per-cpu=2g                 # mem per CPU
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=BLAST_%j.out     # Standard output and error log
#SBATCH --error=BLAST_%j.err		# Standard error log
#SBATCH --qos=page-b

pwd; hostname; date


module load python/2.7.14

python 14a_ortholog_filter.py 1by1.out B_platifrons
