#!/bin/bash

conda activate ldsc

#### gzip all the input bed files first

#gzip /projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/ldsc_multiomeRNA_DLPFC/marker_gene_inputs_to_ldsc/foreground/introns_and_flanks/*
#gzip /projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/ldsc_multiomeRNA_DLPFC/marker_gene_inputs_to_ldsc/background/introns_and_flanks/*

PathToBed="/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/ldsc_multiomeRNA_DLPFC/marker_gene_inputs_to_ldsc/foreground/introns_and_flanks/*.bed.gz"
PathToBackgroundBed="/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/ldsc_multiomeRNA_DLPFC/marker_gene_inputs_to_ldsc/background/introns_and_flanks/*.bed.gz"

EnrichmentsDir='/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/ldsc_multiomeRNA_DLPFC/gwas_enrichments/'

if [ ! -d "$EnrichmentsDir" ]; then
  mkdir "$EnrichmentsDir"
  echo "made $EnrichmentsDir"
fi

BackgroundDir=${EnrichmentsDir}"backgroundBedFiles"
if [ ! -d "$BackgroundDir" ]; then
  mkdir "$BackgroundDir"
  echo "made $BackgroundDir"
fi

PathToOutput=${EnrichmentsDir}"backgroundBedFiles/all_cell_types.bed"

echo "Making union of foregrounds"

### make the union of all foregrounds - should take a few seconds
zcat ${PathBed} | awk 'OFS="\t" {print $1, $2, $3}'| bedtools sort -i stdin | bedtools merge -i stdin > ${PathToOutput}
gzip ${PathToOutput}

echo "Done. Merging union of foregrounds with all genes"

### merge the union of all foregrounds with DHS - should take a few seconds

PathToMerged=${EnrichmentsDir}"backgroundBedFiles/all_cell_types.bed.gz"

# Merge - basically pointless, but it's fast
PathToOutput=${EnrichmentsDir}"backgroundBedFiles/all_genes.bed"
zcat ${PathToMerged} ${PathToBackgroundBed} | awk 'OFS="\t" {print $1, $2, $3}' | bedtools sort -i stdin | bedtools merge -i stdin > ${PathToOutput}


echo "Done. Annotating background"

### annotate the background - should take a few hours (around 3)
Annotations=${EnrichmentsDir}"annotations/"
if [ ! -d "$Annotations" ]; then
  mkdir "$Annotations"
  echo "made $Annotations"
fi

PathToAnnotOutput=${EnrichmentsDir}"annotations/backgroundAnnotations/"
if [ ! -d "$PathToAnnotOutput" ]; then
  mkdir "$PathToAnnotOutput"
  echo "made $PathToAnnotOutput"
fi

sbatch --partition pool1 /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh -i ${PathToOutput} -n all_genes -o ${PathToAnnotOutput}

echo "Done. Converting to bed"
### convert narrowPeak foreground files to bed file format - should take a few seconds

PathToOutput=${EnrichmentsDir}"foregroundBedFiles/"
if [ ! -d "$PathToOutput" ]; then
  mkdir "$PathToOutput"
  echo "made $PathToOutput"
fi

cp ${PathToBed} ${PathToOutput}
gunzip ${PathToOutput}*

#for file in ${PathToNarrowPeak}; do x=${file##*/}; cat ${PathToOutput}${x/.narrowPeak.gz/.bed}; done
#for file in ${PathToBed}; do x=${file##*/}; cat ${file} | awk 'OFS="\t" {print $1, $2, $3}' > ${PathToOutput}${x}; done

#echo "Done. Annotating foregrounds"

#### annotate foregrounds - should take a few hours (around 6+)
PathToAnnotOutput=${EnrichmentsDir}"annotations/foregroundAnnotations/"
if [ ! -d "$PathToAnnotOutput" ]; then
  mkdir "$PathToAnnotOutput"
  echo "made $PathToAnnotOutput"
fi

for file in ${PathToOutput}*.bed; do x=${file##*/}; base_file_without_bed=${x/.bed}; sbatch --partition pool1 --mem=5GB /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh -i $file -n $base_file_without_bed -o ${PathToAnnotOutput}; sleep 2m; done

# munge gwas'
# create ldcts file
# run: sbatch file_name.sb