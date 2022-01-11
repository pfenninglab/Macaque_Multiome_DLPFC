#!/bin/bash
#SBATCH --partition=pool1,pfen1
#SBATCH --mem=23G
#SBATCH --export=ALL
#SBATCH --signal=2
#SBATCH --job-name=mk_arrow
#SBATCH --error=logs/mk_arrow_%A_%a_out.txt
#SBATCH --output=logs/mk_arrow_%A_%a_out.txt
#SBATCH --array 1-5
#SBATCH --no-requeue

ulimit -n 16000; source ~/.bashrc
PROJDIR=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC
CODEDIR=$PROJDIR/code/raw_code/preprocess_multiomeATAC_reads
DATADIR=$PROJDIR/data/raw_data; TMPDIR=/scratch/bnphan; BEDDIR=$DATADIR/bed
cd $CODEDIR
mkdir -p $TMPDIR $BEDDIR

# sample-specific files
NAME=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $5}' ${DATADIR}/tables/SampleSheet_ATAC_DLPFC.csv )
TISSUE=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $3}' ${DATADIR}/tables/SampleSheet_ATAC_DLPFC.csv )
BATCH=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $2}' ${DATADIR}/tables/SampleSheet_ATAC_DLPFC.csv )

#########################################
# make the arrow file for this fragments file
ALIGN_BED=${NAME}_${TISSUE}-${BATCH}.corrected.tsv.gz
LABEL=${NAME}_${TISSUE}-${BATCH}

if [[ ! -f ${BEDDIR}/${ALIGN_BED} ]]; then
	echo "File not found: ${BEDDIR}/${ALIGN_BED}."
elif [[ -f ${ARROWDIR}/${LABEL}.arrow ]]; then
	echo "Arrow file already made: ${ARROWDIR}/${LABEL}.arrow"
else
	echo "Creating arrow file now: ${ARROWDIR}/${LABEL}.arrow"
	source activate r4
	cd $TMPDIR
	rsync -Paq ${BEDDIR}/${ALIGN_BED}* $TMPDIR
	# make the arrow file
	Rscript --vanilla $CODEDIR/make_ArchR_arrows_rheMac10.R -b ${ALIGN_BED} -s ${LABEL}
	rsync --remove-source-files -Paq ${LABEL}.arrow ${ARROWDIR}

	# copy the QC folder and clean up
	rsync --remove-source-files -Paq QualityControl/${LABEL} ${ARROWDIR}/QualityControl
	rsync --remove-source-files -Paq ArchRLogs/ ${ARROWDIR}/ArchRLogs
	rm ${ALIGN_BED}*
fi


