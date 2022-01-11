#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time 3-00:00:00
#SBATCH --mem=15G
#SBATCH --job-name=meme
#SBATCH --error=logs/meme-chip_%A_%a_out.txt
#SBATCH --output=logs/meme-chip_%A_%a_out.txt
#SBATCH --array=1-23%4

PROJDIR=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/raw_data
DBFILE=/home/bnphan/resources/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
GENOME=/home/bnphan/resources/genomes/rheMac10/rheMac10.fa
TMPDIR=/scratch/bnphan/meme-chip; mkdir -p $TMPDIR

########################################
# 1) get the narrowpeak file
NARROWPEAK=$(ls -l ${PROJDIR}/peak_rheMac10/Macaque_PFC* | \
	awk -v IND="${SLURM_ARRAY_TASK_ID}" 'NR==(IND + 1) {print $9}')
TMPBED=${TMPDIR}/$( basename $NARROWPEAK | sed 's/.narrowPeak.gz/.narrowPeak/g' )
zcat $NARROWPEAK > ${TMPBED}

##########################################
# 2) extract fasta sequences for cell type
POSFILE=${TMPDIR}/$( basename $NARROWPEAK | sed 's/.narrowPeak.gz/.fasta/g')
bedtools getfasta -name -fi ${GENOME} -bed ${TMPBED} > ${POSFILE} 

############################################################
# 3) create negative sequences against all peaks in PFC
BGDNARROWPEAK=${PROJDIR}/peak_rheMac10/Macaque_PFC_merged_peaks.bgd.narrowPeak
if [[ ! -f $BGDNARROWPEAK ]]; then zcat ${PROJDIR}/peak_rheMac10/Macaque_PFC* > $BGDNARROWPEAK; fi
NEGFILE=${PROJDIR}/peak_rheMac10/Macaque_PFC_merged_peaks.bgd.fasta
if [[ ! -f $NEGFILE ]]; then bedtools getfasta -name -fi ${GENOME} -bed ${BGDNARROWPEAK} > ${NEGFILE}; fi

######################################
# 4) make output dir and run meme-chip 
OUTDIR=${PROJDIR}/meme/$(echo $NARROWPEAK | sed 's/.narrowPeak.gz/vsMergePeaks/g' )
mkdir -p $OUTDIR
/home/ikaplow/anaconda2/bin/meme-chip \
-fimo-skip -spamo-skip -noecho \
-oc $OUTDIR -db $DBFILE -neg $NEGFILE $POSFILE

rm ${TMPBED} ${POSFILE} ${OUTDIR}/*.fasta