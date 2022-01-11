#################
## directories ##
SETDIR=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
CODEDIR=${SETDIR}/code/raw_code/ldsc_multiomeATAC_DLPFC
DATADIR=${SETDIR}/data/tidy_data/ldsc_multiomeATAC_DLPFC
ANNOTDIR=${DATADIR}/annotations
cd $CODEDIR; mkdir -p logs $ANNOTDIR; mv -n annot* logs

#######################################
## merge human epigenome backgrounds ##
BG1=${SETDIR}/data/raw_data/halper
BG2=/home/bnphan/src/atac-seq-pipeline/genome/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz
BGNAME=Macaque_PFC.BG_Honeybadger
BGFILE=${DATADIR}/${BGNAME}.bed.gz

## make the background peak file
if [[ ! -f ${BGFILE} ]]; then
cat ${BG1}/Macaque_PFC*Macaca_mulattaToHomo_sapiens.HALPER.narrowPeak.gz ${BG2} | zcat | \
cut -f 1-3 | sort --parallel=10 -k1,1 -k2,2n | gzip > ${DATADIR}/${BGNAME}.tmp.bed.gz 
bedtools merge -i ${DATADIR}/${BGNAME}.tmp.bed.gz | gzip > ${BGFILE}
rm ${DATADIR}/${BGNAME}.tmp.bed.gz
fi

#######################################
# for background for binary annotations
if [[ ! -f "${ANNOTDIR}/${BGNAME}.AFR.1.l2.ldscore.gz" ]]; then 
	sbatch --mem 2G -p pfen1,pfen_bigmem,pfen3,pfen2 --job-name=${BGNAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
	-i ${BGFILE} -n ${BGNAME} -g hg38 -p AFR -o $ANNOTDIR
fi
if [[ ! -f "${ANNOTDIR}/${BGNAME}.EUR.1.l2.ldscore.gz" ]]; then 
	sbatch --mem 2G -p pfen1,pfen_bigmem,pfen3,pfen2 --job-name=${BGNAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
	-i ${BGFILE} -n ${BGNAME} -g hg38 -p EUR -o $ANNOTDIR
fi

#############################################
## annotate for LDSC w/ binary annotations ##
CTS_AFR_FN=${DATADIR}/Macaque_PFC.AFR_hg38.ldcts; > $CTS_AFR_FN
CTS_EUR_FN=${DATADIR}/Macaque_PFC.EUR_hg38.ldcts; > $CTS_EUR_FN
for BED in ${BG1}/*Macaca_mulattaToHomo_sapiens.HALPER.narrowPeak.gz ; do
NAME=$(basename $BED | sed 's/.GenBankRheMac8.Macaca_mulattaToHomo_sapiens.HALPER.narrowPeak.gz//g')
if [[ ! -f "${ANNOTDIR}/${NAME}.AFR.1.l2.ldscore.gz" ]]; then 
echo "Annotations for ${NAME} not found."
sbatch --mem 2G -p pfen1,pfen_bigmem,pfen3,pfen2 --job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
fi
if [[ ! -f "${ANNOTDIR}/${NAME}.EUR.1.l2.ldscore.gz" ]]; then 
echo "Annotations for ${NAME} not found."
sbatch --mem 2G -p pfen1,pfen_bigmem,pfen3,pfen2 --job-name=${NAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
fi
echo -e "${NAME}\t${ANNOTDIR}/${NAME}.AFR.,${ANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
echo -e "${NAME}\t${ANNOTDIR}/${NAME}.EUR.,${ANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN
done


