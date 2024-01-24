#!/bin/bash
#SBATCH --job-name=igDetect    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20                    # Run on a single CPU 
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/igDetective%j.log   # Standard output and error log
#SBATCH --mem=60G




source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/IGdetective
conda env list

/home1/zhuyixin/.conda/envs/IGdetective/bin/python /home1/zhuyixin/IgDetective/run_iterative_igdetective.py $1 $2
