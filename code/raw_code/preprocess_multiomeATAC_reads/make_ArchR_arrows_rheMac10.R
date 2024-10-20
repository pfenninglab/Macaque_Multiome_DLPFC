suppressMessages(library(ArchR))

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

#####################################
### read in command line options ####
library(optparse)
option_list = list(
  make_option(c("-b", "--FileName"), type="character", default=NULL, 
              help="bam/bed/fragments file to make into arrowfile", metavar="character"),
  make_option(c("-s", "--SampleName"), type="character", default=NULL, 
              help="arrow file sample name", metavar="character"),
  make_option(c("-f", "--force"), type="logical", default=FALSE, 
              help="whether to force making this file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$FileName)){
  print_help(opt_parser)
  stop("Please provide Fragments File name, SampleName, and outArrowName.n", call.=FALSE)
}

if (is.null(opt$SampleName)){
  # use base file to name arrows if not provided
  opt$SampleName = gsub('.bam|.bed.gz|tsv.gz',basename(opt$SampleName))
}

########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = round(parallel::detectCores()/2))
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
# contains `geneAnnotation` and `genomeAnnotation` objects
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))

##############################
### make that arrow file! ####
ArrowFiles <- createArrowFiles( inputFiles = opt$FileName, 
                                sampleNames = opt$SampleName, 
                                minTSS = 4, #Dont set this too high because you can always increase later
                                minFrags = 1000, addTileMat = TRUE,
                                addGeneScoreMat = TRUE,
                                geneAnnotation = geneAnnotation, #must be the custom rheMac10 version
                                genomeAnnotation = genomeAnnotation, #must be the custom rheMac10 version
                                cleanTmp = FALSE,
                                force = FALSE)

doubleScores <- addDoubletScores(input = ArrowFiles, 
                                k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                                LSIMethod = 1,  knnMethod = "UMAP") #Refers to the embedding to use for nearest neighbor search with doublet projection.

