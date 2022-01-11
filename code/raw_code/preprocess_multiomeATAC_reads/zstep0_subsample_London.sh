## the data transfer through BaDoi's google drive

## ATAC for monkey London
cd /projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/raw_data/fastq/ATAC/London

mkdir -p ~/subsample_fastq/ATAC

for FILE in *_L00[1-2]_R?_001.fastq.gz; do 
FILE2=~/subsample_fastq/ATAC/$(basename $FILE .fastq.gz).subSample1M.fastq.gz
zcat $FILE | head -1000000 | gzip > $FILE2
done

cd /projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/raw_data/fastq/RNA/London

mkdir -p ~/subsample_fastq/RNA

for FILE in *_L00[1-2]_??_001.fastq.gz; do 
FILE2=~/subsample_fastq/RNA/$(basename $FILE .fastq.gz).subSample1M.fastq.gz
zcat $FILE | head -1000000 | gzip > $FILE2
done
