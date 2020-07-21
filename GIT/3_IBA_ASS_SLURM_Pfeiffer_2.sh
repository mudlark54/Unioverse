#!/bin/sh
#SBATCH --job-name=IBA_hog	 	# Job name 			
#SBATCH --mail-type=ALL             	 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=pfeifferj@si.edu 	 # Where to send mail	
#SBATCH --nodes=1 	
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem-per-cpu=2g                 # mem per CPU
#SBATCH --time=60:00:00              # Time limit hrs:min:sec
#SBATCH --output=IBA_hog_%j.out     # Standard output and error log
#SBATCH --error=IBA_hog_%j.err		# Standard error log
#SBATCH --qos=page-b
#SBATCH --array=1-26

INFILE1=$(head -n $SLURM_ARRAY_TASK_ID RAWLIST | tail -n1|cut -f1 | sed "s/.fastq.gz/_val_1.fq/")
INFILE2=$(head -n $SLURM_ARRAY_TASK_ID RAWLIST | tail -n1|cut -f2 | sed "s/.fastq.gz/_val_2.fq/")
TAXA=$(head -n $SLURM_ARRAY_TASK_ID RAWLIST| tail -n1|cut -f3)
#echo "$SLURM_ARRAY_TASK_ID $TAXA" >> genelist


mkdir ./IBAass_$TAXA
cd ./IBAass_$TAXA

module load usearch/11.0.667
module load python/3.8
module load bridger/20141201

export OMP_NUM_THREADS=2
python ../IBA.py -raw1 ../$INFILE1 -raw2 ../$INFILE2 -d ../REF -n 2 -t 1 -p 16 -g 200 -c 2 -label B_platifrons -k 25 -taxa $TAXA


