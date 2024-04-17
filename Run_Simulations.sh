#!/bin/bash

#SBATCH -o Protein9_Result/Replicate_%a.Rout
#SBATCH --array=1-2000
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --mem-per-cpu=3gb
#SBATCH -t 48:00:00

module load R

R CMD BATCH --vanilla Protein9_Main.R  Protein9_Result/Replicate_${SLURM_ARRAY_TASK_ID}.Rout
