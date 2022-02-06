#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --time 3-00:00:00
#SBATCH --mem=23G
#SBATCH --job-name=meme
#SBATCH --error=logs/meme-diff_%A_%a_out.txt
#SBATCH --output=logs/meme-diff_%A_%a_out.txt
#SBATCH --array=1-91%20

## commenting for mike
PROJDIR=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC
CODEDIR=$PROJDIR/code/raw_code/celltype_specific_enhancers
DATADIR=$PROJDIR/data/tidy_data/celltype_specific_enhancers
DBFILE=/home/bnphan/resources/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
GENOME=/home/bnphan/resources/genomes/rheMac10/rheMac10.fa
cd $CODEDIR

##################################
### cell type and fold index #####
MODEL=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = "," } ;NR==IND + 1 {print $1}' $DATADIR/tables/cell_type_models_to_train_DLPFC.csv`
CELLTYPE=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = "," } ;NR==IND + 1 {print $2}' $DATADIR/tables/cell_type_models_to_train_DLPFC.csv`
TRAINPOSFILE=$DATADIR/fasta/rheMac10_${MODEL}_fold1_train_positive.fa
TRAINNEGFILE=$DATADIR/fasta/rheMac10_${MODEL}_fold1_train_negative.fa

######################################
# 4) make output dir and run meme-chip 
OUTDIR=${DATADIR}/meme/${CELLTYPE}/${MODEL}
mkdir -p $OUTDIR

if [[ ! -f ${OUTDIR}/meme-chip.html ]]; then
/home/ikaplow/anaconda2/bin/meme-chip -fimo-skip -spamo-skip -noecho  \
-oc $OUTDIR -db $DBFILE -neg $TRAINNEGFILE $TRAINPOSFILE
rm ${OUTDIR}/*.fa  ${OUTDIR}/seqs* ${OUTDIR}/control-centered
else
echo "already completed meme-chip run for ${OUTDIR}."
fi
