#!/bin/sh
#SBATCH --job-name=userach	 	# Job name 			
#SBATCH --mail-type=ALL             	 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jpfeiffer@ufl.edu 	 # Where to send mail	
#SBATCH --nodes=1 	
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem-per-cpu=2g                 # mem per CPU
#SBATCH --time=60:00:00              # Time limit hrs:min:sec
#SBATCH --output=IBA_hog_%j.out     # Standard output and error log
#SBATCH --error=IBA_hog_%j.err		# Standard error log
#SBATCH --qos=page-b


module load usearch/11.0.667

usearch -fastx_getseqs ALL_FULL_LOCI_w_ref_sl.fa -labels output1_keep.list -label_substr_match -fastaout FINAL_FULL_CLEAN.fa
