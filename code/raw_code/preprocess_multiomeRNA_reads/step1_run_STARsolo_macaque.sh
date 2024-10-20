#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen3
#SBATCH --job-name=STARsolo
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --error=logs/STARsolo_%A_%a.txt
#SBATCH --output=logs/STARsolo_%A_%a.txt
#SBATCH --array=1-2

source ~/.bashrc

###########################################
# get sample name in demultiplex csv table
PROJDIR=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC
DATADIR=${PROJDIR}/data/raw_data
FASTQDIR=$DATADIR/fastq/RNA
SOLODIR=$DATADIR/STARsolo_out

mkdir -p $TMPDIR $SOLODIR; cd $SOLODIR

# sample-specific files
DIRLABEL=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $1}' ${DATADIR}/tables/SampleSheet_RNA.csv )
SAMPLE=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $4}' ${DATADIR}/tables/SampleSheet_RNA.csv )
PREFIX=${DIRLABEL}_DLPFC-1.

##################################################
# cDNA fragment in Read2, cell barcode in Read1
cDNA_FASTQ=$(ls ${FASTQDIR}/${DIRLABEL}/*R2_001.fastq.gz | tr '[[:space:]]' ',' | sed 's/,$//g')
CB_FASTQ=$(ls ${FASTQDIR}/${DIRLABEL}/*R1_001.fastq.gz | tr '[[:space:]]' ',' | sed 's/,$//g')

# perform mapping & read quantification
GENOME_DIR=/home/bnphan/resources/genomes/rheMac10
BARCODES=/home/bnphan/resources/cell_ranger_barcodes/737K-arc-v1_RNA.txt

if [[ ! -f "${PREFIX}Log.final.out" ]]; then
	echo "Aligning samples w/ STARsolo for: ${SAMPLE}."
	~/src/STAR-2.7.9a/bin/Linux_x86_64/STAR \
	--outFileNamePrefix $PREFIX \
	--soloType Droplet \
	--readFilesIn $cDNA_FASTQ $CB_FASTQ \
	--readFilesCommand zcat \
	--genomeDir $GENOME_DIR \
	--limitOutSJcollapsed 5000000 \
	--runThreadN 12 \
	--soloCBwhitelist $BARCODES \
	--soloFeatures GeneFull \
	--soloBarcodeReadLength 0 \
	--soloCBmatchWLtype 1MM \
	--soloCellFilter EmptyDrops_CR \
	--soloMultiMappers EM \
	--soloUMIlen 12 \
	--soloUMIdedup 1MM_CR \
	--outSAMtype None
	# --outSAMtype BAM SortedByCoordinate \
	# --outSJtype None \
	# --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM CY UY
	rm -rf ${PREFIX}_STARtmp
else 
	echo "A completed STARsolo already exists: ${SOLODIR}/${PREFIX}Solo.out"
fi
