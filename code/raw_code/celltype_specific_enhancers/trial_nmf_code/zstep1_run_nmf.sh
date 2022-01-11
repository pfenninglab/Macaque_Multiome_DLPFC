#!/bin/bash
#SBATCH --partition=pool1
#SBATCH --mem=23G
#SBATCH --time 1-00:00:00
#SBATCH --export=ALL
#SBATCH --signal=2
#SBATCH --job-name=run_nmf
#SBATCH --error=logs/run_nmf_%A_%a_out.txt
#SBATCH --output=logs/run_nmf_%A_%a_out.txt
#SBATCH --array 1-100%20
#SBATCH --no-requeue

ulimit -n 16000; source ~/.bashrc; source activate r4
PROJDIR=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC
CODEDIR=$PROJDIR/code/raw_code/celltype_specific_enhancers
cd $CODEDIR

Rscript --vanilla $CODEDIR/NMF_factor_peakMatrix_run.R -i ${SLURM_ARRAY_TASK_ID}