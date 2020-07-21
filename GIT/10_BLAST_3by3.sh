#!/bin/sh
#SBATCH --job-name=BLAST	 	# Job name 			
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


module load ncbi_blast/2.6.0


blastn -task blastn -query 3by3.in -db B_platifrons -out 3by3.out -outfmt 6 -max_target_seqs 3 -max_hsps 3 -num_threads 8

