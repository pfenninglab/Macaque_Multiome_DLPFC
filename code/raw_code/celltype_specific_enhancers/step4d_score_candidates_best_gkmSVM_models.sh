#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1,pool3-bigmem,pfen1
#SBATCH --time=3-00:00:00
#SBATCH --mem=8G
#SBATCH --job-name=svm_pfc
#SBATCH --error=logs/lsgkm_%A_%a_output.txt
#SBATCH --output=logs/lsgkm_%A_%a_output.txt
#SBATCH --array=1-87%20

PROJDIR=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC
DATADIR=$PROJDIR/data/tidy_data/celltype_specific_enhancers
cd $PROJDIR/code/raw_code/celltype_specific_enhancers

######################################################
# input fasta files for training, validation, and test sets
TESTPOSFILE=$DATADIR/fasta/rheMac10_${LABEL}_test_positive.fa
TESTNEGFILE=$DATADIR/fasta/rheMac10_${LABEL}_test_negative.fa
CANDIDATES=$DATADIR/fasta/rheMac10_DLPFC.noncoding_peak.fa

###########################
# cell type and fold index
TRAINED_MODELS=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = "\t" } ;NR==IND + 1 {print $12}' $DATADIR/tables/cell_type_bestSVMs_DLPFC.tsv`
CELLTYPE=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = "\t" } ;NR==IND + 1 {print $2}' $DATADIR/tables/cell_type_bestSVMs_DLPFC.tsv`
PREDICTBASE=$DATADIR/predictions/${CELLTYPE}/$( basename $TRAINED_MODELS .model.txt )

#####################################################################
# score the set of candidate enhancer sequences w/ the best models ##
FILE_C=${PREDICTBASE}_gkmpredict_noncoding_peak.txt
if [[ ! -f $FILE_C ]]; then gkmpredict -T 4 $CANDIDATES ${TRAINED_MODELS} ${FILE_C}; fi
# # FILE=${PREDICTBASE}_gkmpredict_test_positive.txt
# # if [[ ! -f $FILE || ! -s $FILE ]]; then gkmpredict -T 4 $TESTPOSFILE ${TRAINED_MODELS} ${PREDICTBASE}_gkmpredict_test_positive.txt; fi
# # FILE=${PREDICTBASE}_gkmpredict_test_positive.txt
# # if [[ ! -f $FILE || ! -s $FILE ]]; then gkmpredict -T 4 $TESTNEGFILE ${TRAINED_MODELS} ${PREDICTBASE}_gkmpredict_test_negative.txt; fi
