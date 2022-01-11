#!/bin/bash
#SBATCH --partition=pfen1
#SBATCH --mem=24G
#SBATCH --export=ALL
#SBATCH --signal=2
#SBATCH --job-name=correct_multiome_barcodes
#SBATCH --error=logs/correct_multiome_barcodes_%A_%a_out.txt
#SBATCH --output=logs/correct_multiome_barcodes_%A_%a_out.txt
#SBATCH --array 1-2
#SBATCH --no-requeue

ulimit -n 16000; source ~/.bashrc
PROJDIR=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC
CODEDIR=$PROJDIR/code/raw_code/preprocess_multiomeATAC_reads
DATADIR=$PROJDIR/data/raw_data; TMPDIR=/scratch/bnphan
BEDDIR=$DATADIR/bed
mkdir -p $TMPDIR $BEDDIR
cd $CODEDIR

# sample-specific files
DIRLABEL=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $1}' ${DATADIR}/tables/SampleSheet.csv )
SAMPLE=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $4}' ${DATADIR}/tables/SampleSheet.csv )

# rheMac10 reference genome indexed for chromap
BARCODE_LIST=/home/bnphan/resources/cell_ranger_barcodes/737K-arc-v1_ATAC2RNA.txt
GENOME_FASTA=/home/bnphan/resources/genomes/rheMac10/rheMac10.fa
CHROMAP_IDX=/home/bnphan/resources/genomes/chromap/rheMac10/index

ALIGN_BED=${DIRLABEL}_${SAMPLE}.aln.bed.gz
ALIGN_BED2=${DIRLABEL}_${SAMPLE}.corrected.tsv

if [[ ! -f ${BEDDIR}/${ALIGN_BED2}.gz && ! -f ${BEDDIR}/${ALIGN_BED2}.gz.tbi ]]; then
cd $BEDDIR
echo "Converting Multiome ATAC cell barcode to Multiome RNA cell barcode for:${DIRLABEL}_${SAMPLE}."
awk 'BEGIN { FS = OFS = "\t"; } 
FNR==NR{a[$1]=$2;next} 
{if ($4 in a)
	print $1,$2,$3,a[$4],$5
else 
	print $0
}' $BARCODE_LIST <(zcat ${ALIGN_BED}) > $ALIGN_BED2
bgzip -c ${ALIGN_BED2} > ${ALIGN_BED2}.gz
tabix -p bed ${ALIGN_BED2}.gz
rm ${ALIGN_BED2}
else 
echo 'Multiome barcode-corrected fragments-file already exists.'
fi



