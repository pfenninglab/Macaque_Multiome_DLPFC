#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1,pool3-bigmem,pfen1
#SBATCH --time=3-00:00:00
#SBATCH --mem=23G
#SBATCH --job-name=svm_pfc
#SBATCH --error=logs/lsgkm_%A_%a_output.txt
#SBATCH --output=logs/lsgkm_%A_%a_output.txt
#SBATCH --array=1-91%40

PROJDIR=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC
DATADIR=$PROJDIR/data/tidy_data/celltype_specific_enhancers
cd $PROJDIR/code/raw_code/celltype_specific_enhancers

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
TRAINPOSFILE=$DATADIR/fasta/rheMac10_${LABEL}_train_positive.fa
TRAINNEGFILE=$DATADIR/fasta/rheMac10_${LABEL}_train_negative.fa
VALIDPOSFILE=$DATADIR/fasta/rheMac10_${LABEL}_valid_positive.fa
VALIDNEGFILE=$DATADIR/fasta/rheMac10_${LABEL}_valid_negative.fa
TESTPOSFILE=$DATADIR/fasta/rheMac10_${LABEL}_test_positive.fa
TESTNEGFILE=$DATADIR/fasta/rheMac10_${LABEL}_test_negative.fa
CANDIDATES=$DATADIR/fasta/rheMac10_DLPFC.noncoding_peak.fa

###############################################
########### train the SVM models ##############
for CPARAM in $(seq .1 0.1 1); do 
echo "Using regulariation parameter ${CPARAM}.";
# train SVM model w/ regular WPARAM and CPARAM
OUTPREFIX=$DATADIR/models_svm/${CELLTYPE}/${LABEL}_t4_l7_k6_d1_c${CPARAM}_w${WPARAM}
if [[ ! -f "${OUTPREFIX}.model.txt" || $TRAINPOSFILE -nt "${OUTPREFIX}.model.txt" ]]; then
gkmtrain -t 4 -l 7 -k 6 -d 1 -c ${CPARAM} -w ${WPARAM} \
-M 50 -H 50 -m 10000 -s -T 4 $TRAINPOSFILE $TRAINNEGFILE $OUTPREFIX
else
echo "Previously trained model found: ${OUTPREFIX}.model.txt"
fi

# train a model w/ the inverse CPARAM for the very imbalanced models
CPARAM2=$(echo "1 / ${CPARAM}"|bc)
OUTPREFIX=$DATADIR/models_svm/${CELLTYPE}/${LABEL}_t4_l7_k6_d1_c${CPARAM2}_w${WPARAM}
if [[ ! -f "${OUTPREFIX}.model.txt" || $TRAINPOSFILE -nt "${OUTPREFIX}.model.txt" ]]; then
gkmtrain -t 4 -l 7 -k 6 -d 1 -c ${CPARAM2} -w ${WPARAM} \
-M 50 -H 50 -m 10000 -s -T 4 $TRAINPOSFILE $TRAINNEGFILE $OUTPREFIX
else
echo "Previously trained model found: ${OUTPREFIX}.model.txt"
fi
done

#############################################################################
# pedict validation and candidate enhancer sequences on same set of models ##
for TRAINED_MODELS in $(ls $DATADIR/models_svm/${CELLTYPE}/${LABEL}*); do 
PREDICTBASE=$DATADIR/predictions/${CELLTYPE}/$( basename $TRAINED_MODELS .model.txt )
FILE_P=${PREDICTBASE}_gkmpredict_valid_positive.txt
if [[ ! -f $FILE_P ]]; then gkmpredict -T 4 $VALIDPOSFILE ${TRAINED_MODELS} ${PREDICTBASE}_gkmpredict_valid_positive.txt; fi
FILE_N=${PREDICTBASE}_gkmpredict_valid_negative.txt
if [[ ! -f $FILE_N ]]; then gkmpredict -T 4 $VALIDNEGFILE ${TRAINED_MODELS} ${PREDICTBASE}_gkmpredict_valid_negative.txt; fi

## evaluate the validation set predictions
SCORE_FN=${PREDICTBASE}_gkmpredict_valid_eval.feather
if [[ -f $FILE_P && -f $FILE_N && ! -f $SCORE_FN ]]; then
python score_gkm_classifier.py --prefix=$CELLTYPE --label=valid --out_dir=$DATADIR \
--model_name=$TRAINED_MODELS --pred_pos=$FILE_P --pred_neg=$FILE_N
fi
# # FILE=${PREDICTBASE}_gkmpredict_test_positive.txt
# # if [[ ! -f $FILE || ! -s $FILE ]]; then gkmpredict -T 4 $TESTPOSFILE ${TRAINED_MODELS} ${PREDICTBASE}_gkmpredict_test_positive.txt; fi
# # FILE=${PREDICTBASE}_gkmpredict_test_positive.txt
# # if [[ ! -f $FILE || ! -s $FILE ]]; then gkmpredict -T 4 $TESTNEGFILE ${TRAINED_MODELS} ${PREDICTBASE}_gkmpredict_test_negative.txt; fi
done
