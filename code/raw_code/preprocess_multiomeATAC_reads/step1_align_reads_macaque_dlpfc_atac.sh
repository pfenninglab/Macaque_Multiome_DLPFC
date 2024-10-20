#!/bin/bash
#SBATCH --partition=pool1,pfen1
#SBATCH --time=2-00:00:00
#SBATCH --mem=23G
#SBATCH --export=ALL
#SBATCH --signal=2
#SBATCH --job-name=chromap
#SBATCH --error=logs/chromap_%A_%a_out.txt
#SBATCH --output=logs/chromap_%A_%a_out.txt
#SBATCH --array 1-5
#SBATCH --no-requeue

ulimit -n 16000; source ~/.bashrc
PROJDIR=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC
CODEDIR=$PROJDIR/code/raw_code/preprocess_multiomeATAC_reads
DATADIR=$PROJDIR/data/raw_data; TMPDIR=/scratch/bnphan; BEDDIR=$DATADIR/bed
cd $CODEDIR
mkdir -p $TMPDIR $BEDDIR

# sample-specific files
SAMPLE=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $6}' ${DATADIR}/tables/SampleSheet_ATAC_DLPFC.csv )
NAME=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $5}' ${DATADIR}/tables/SampleSheet_ATAC_DLPFC.csv )
TISSUE=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $3}' ${DATADIR}/tables/SampleSheet_ATAC_DLPFC.csv )
BATCH=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $2}' ${DATADIR}/tables/SampleSheet_ATAC_DLPFC.csv )
CHEM=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $9}' ${DATADIR}/tables/SampleSheet_ATAC_DLPFC.csv )

# rheMac10 reference genome indexed for chromap
if [[ $CHEM == "Multiome kit" ]]; then
FASTQDIR=$DATADIR/fastq/ATAC
BARCODE_LIST=/home/bnphan/resources/cell_ranger_barcodes/737K-arc-v1_ATAC_revComp.txt
else # regular snATAC-seq barcode list
BARCODE_LIST=/home/bnphan/resources/cell_ranger_barcodes/737K-cratac-v1.txt
FASTQDIR=/data/pfenninggroup/singleCell/Macaque_snATAC-seq/m015_Peanut/snATACseq
fi
GENOME_FASTA=/home/bnphan/resources/genomes/rheMac10/rheMac10.fa
CHROMAP_IDX=/home/bnphan/resources/genomes/chromap/rheMac10/index
ALIGN_BED=${NAME}_${TISSUE}-${BATCH}.aln.bed

##################################
# align scATAC reads w/ chromap
if [[ ! -f ${BEDDIR}/${ALIGN_BED}.gz && ! -f ${BEDDIR}/${ALIGN_BED}.gz.tbi ]]; then
cd $TMPDIR; 

## get the fastq files
echo "1. Copying fastqs for: ${NAME}_${TISSUE}-${BATCH}."
rsync -Paq ${FASTQDIR}/*/${SAMPLE}*_R?_001.fastq.gz .

FQ1=$(ls ${SAMPLE}*_R1_001.fastq.gz | tr '\n' ',' | sed 's/,$//g')
FQ2=$(ls ${SAMPLE}*_R3_001.fastq.gz | tr '\n' ',' | sed 's/,$//g' )

## make the BC file the rev complement and trim spacer for multiome files
if [[ $CHEM == "Multiome kit" ]]; then
echo "2. Detected multiomeATAC sample."
echo "Creating the trimmed ARC cell barcode for: ${NAME}_${TISSUE}-${BATCH}."
for BC in $(ls ${SAMPLE}*_R2_001.fastq.gz); do
BCX=$(basename $BC _R2_001.fastq.gz)_R2_trimmed.fastq.gz
if [ ! -f $BCX ]; then seqtk trimfq -b 8 $BC | gzip > $BCX; fi
done
BC1=$(ls ${SAMPLE}*_R2_trimmed.fastq.gz | tr '\n' ',' | sed 's/,$//g')
else # regular snATAC-seq samples, don't need barcode trimming
echo "2. Detected snATAC sample: ${NAME}_${TISSUE}-${BATCH}."
BC1=$(ls ${SAMPLE}*_R2_001.fastq.gz | tr '\n' ',' | sed 's/,$//g')
fi

## align w/ chromap
echo "3. Aligning scATAC samples w/ chromap for: ${NAME}_${TISSUE}-${BATCH}."
chromap --preset atac -t 8 \
-x $CHROMAP_IDX -r ${GENOME_FASTA} \
--remove-pcr-duplicates-at-cell-level \
--bc-error-threshold 2 \
--output-mappings-not-in-whitelist \
--bc-probability-threshold .9 \
-1 ${FQ1} -2 ${FQ2} -b ${BC1} -o ${ALIGN_BED} \
--barcode-whitelist ${BARCODE_LIST}

## move trimmed fastq to de-multiplex dir
bgzip -@ 8 -c ${ALIGN_BED} > ${ALIGN_BED}.gz
tabix -p bed ${ALIGN_BED}.gz
rsync -Paq --remove-source-files ${ALIGN_BED}.gz* ${BEDDIR}
rm -f ${SAMPLE}*.fastq.gz ${ALIGN_BED} # clean up any old trim files

else echo 'Raw chromap fragments file already exists.'
fi

#####################################################
## Correct Multiome ATAC cell barcodes to RNA barcode
ALIGN_BED2=${NAME}_${TISSUE}-${BATCH}.corrected.tsv
BARCODE_LIST2=/home/bnphan/resources/cell_ranger_barcodes/737K-arc-v1_ATAC2RNA.txt

if [[ -f ${BEDDIR}/${ALIGN_BED2}.gz && ! -f ${BEDDIR}/${ALIGN_BED2}.gz.tbi ]]; then
echo 'Barcode-corrected fragments file already exists.'
elif [[ $CHEM == "Multiome kit" ]]; then
cd $BEDDIR
echo "Converting Multiome ATAC cell barcode to Multiome RNA cell barcode for: ${NAME}_${TISSUE}-${BATCH}."
awk 'BEGIN { FS = OFS = "\t"; } FNR==NR{a[$1]=$2;next} {if ($4 in a) {print $1,$2,$3,a[$4],$5} else {print $0}}' $BARCODE_LIST2 <(zcat ${ALIGN_BED}.gz ) > $ALIGN_BED2
bgzip -@ 8 -c ${ALIGN_BED2} > ${ALIGN_BED2}.gz
tabix -p bed ${ALIGN_BED2}.gz
rm ${ALIGN_BED2} 
else # if regular old ATAC-seq files
cp ${ALIGN_BED}.gz ${ALIGN_BED2}.gz
cp ${ALIGN_BED}.gz.tbi ${ALIGN_BED2}.gz.tbi 
fi

## clean up unused fragments file
rm ${ALIGN_BED}.gz ${ALIGN_BED}.gz.tbi