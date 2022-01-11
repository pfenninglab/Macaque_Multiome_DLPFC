#!/bin/bash
#SBATCH --partition=pfen1
#SBATCH --mem=24G
#SBATCH --export=ALL
#SBATCH --signal=2
#SBATCH --job-name=chromap
#SBATCH --error=logs/chromap_%A_%a_out.txt
#SBATCH --output=logs/chromap_%A_%a_out.txt
##SBATCH --array 1-2
#SBATCH --no-requeue

ulimit -n 16000; source ~/.bashrc
PROJDIR=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC
CODEDIR=$PROJDIR/code/raw_code/preprocess_scATAC_reads
DATADIR=$PROJDIR/data/raw_data; TMPDIR=/scratch/bnphan
FASTQDIR=$DATADIR/fastq/ATAC
BEDDIR=$DATADIR/bed
mkdir -p $TMPDIR $BEDDIR
cd $CODEDIR

# sample-specific files
DIRLABEL='Fiona'
SAMPLE='1_BYR86A1'

# rheMac10 reference genome indexed for chromap
BARCODE_LIST=/home/bnphan/src/cellranger-atac-1.2.0/cellranger-atac-cs/1.2.0/lib/python/barcodes/Combined_forAndRevComp_737K-cratac-v1.txt
GENOME_FASTA=/home/bnphan/resources/genomes/rheMac10/rheMac10.fa
CHROMAP_IDX=/home/bnphan/resources/genomes/chromap/rheMac10/index
ALIGN_BED=${DIRLABEL}.${SAMPLE}.aln.bed

if [ ! -f ${BEDDIR}/$ALIGN_BED ]; then
	cd $TMPDIR; 

	## get the fastq files
	echo "1. Copying fastqs for: ${SAMPLE}."
	rsync -Paq ${FASTQDIR}/${DIRLABEL}/${SAMPLE}*_??_001.fastq.gz .

	## align w/ chromap
	FQ1=$(ls ${SAMPLE}*_R1_001.fastq.gz | tr '\n' ',' | sed 's/,$//g')
	BC1=$(ls ${SAMPLE}*_R2_001.fastq.gz | tr '\n' ',' | sed 's/,$//g')
	FQ2=$(ls ${SAMPLE}*_R3_001.fastq.gz | tr '\n' ',' | sed 's/,$//g' )
	IDX=$(ls ${SAMPLE}*_I1_001.fastq.gz | tr '\n' ',' | sed 's/,$//g' )

	echo "3. Aligning scATAC samples w/ chromap for: ${SAMPLE}."
	~/src/chromap-0.1_x64-linux/chromap --preset atac \
	-x $CHROMAP_IDX -r ${GENOME_FASTA} -t 8 \
	--remove-pcr-duplicates-at-cell-level \
	--bc-error-threshold 2 \
	--output-mappings-not-in-whitelist \
	--bc-probability-threshold .9 \
	-1 ${FQ1} -2 ${FQ2} -b ${BC1} -o ${ALIGN_BED} \
	--barcode-whitelist ${BARCODE_LIST}

	## move trimmed fastq to de-multiplex dir
	bgzip -c ${ALIGN_BED} > ${ALIGN_BED}.gz
	tabix -p bed ${ALIGN_BED}.gz
	rsync -Paq --remove-source-files ${ALIGN_BED}.gz* ${BEDDIR}
	rm -f ${SAMPLE}*_??_001.fastq.gz ${ALIGN_BED} # clean up any old trim files

else echo 'chromap alignment bed file exists.'
fi



