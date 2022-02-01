#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1,pool3-bigmem,pfen1
#SBATCH --time=1-00:00:00
#SBATCH --mem=1G
#SBATCH --job-name=eval_pfc
#SBATCH --error=logs/score_svm_%A_%a_output.txt
#SBATCH --output=logs/score_svm_%A_%a_output.txt
#SBATCH --array=1-91%40

PROJDIR=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC
DATADIR=$PROJDIR/data/tidy_data/celltype_specific_enhancers
CODEDIR=$PROJDIR/code/raw_code/celltype_specific_enhancers
cd $CODEDIR

###########################
# cell type and fold index
MODEL=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = "," } ;NR==IND + 1 {print $1}' $DATADIR/tables/cell_type_models_to_train_DLPFC.csv`
CELLTYPE=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = "," } ;NR==IND + 1 {print $2}' $DATADIR/tables/cell_type_models_to_train_DLPFC.csv`
mkdir -p $DATADIR/models_svm/${CELLTYPE}/
mkdir -p $DATADIR/predictions/${CELLTYPE}/

# weight of the positives to negatives
WPARAM=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = "," } ;NR==IND + 1 {print $10}' $DATADIR/tables/cell_type_models_to_train_DLPFC.csv`
LABEL=${MODEL}_fold1

# input fasta files for training, validation, and test sets
VALIDPOSFILE=$DATADIR/fasta/rheMac10_${LABEL}_valid_positive.fa
VALIDNEGFILE=$DATADIR/fasta/rheMac10_${LABEL}_valid_negative.fa
TESTPOSFILE=$DATADIR/fasta/rheMac10_${LABEL}_test_positive.fa
TESTNEGFILE=$DATADIR/fasta/rheMac10_${LABEL}_test_negative.fa

#############################################################################
# pedict validation and candidate enhancer sequences on same set of models ##
for TRAINED_MODELS in $(ls $DATADIR/models_svm/${CELLTYPE}/${LABEL}*); do 
PREDICTBASE=$DATADIR/predictions/${CELLTYPE}/$( basename $TRAINED_MODELS .model.txt )
FILE_P=${PREDICTBASE}_gkmpredict_valid_positive.txt
FILE_N=${PREDICTBASE}_gkmpredict_valid_negative.txt
SCORE_FN=${PREDICTBASE}_gkmpredict_valid_eval.feather

if [[ -f $FILE_P && -f $FILE_N ]]; then
python score_gkm_classifier.py --prefix=$CELLTYPE --label=valid --out_dir=$DATADIR \
--model_name=$TRAINED_MODELS --pred_pos=$FILE_P --pred_neg=$FILE_N --force
fi
done
