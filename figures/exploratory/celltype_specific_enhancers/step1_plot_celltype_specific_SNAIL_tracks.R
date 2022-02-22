ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(tidyverse))
library(here)

addArchRThreads(threads = 4)

######################
## load ArchR projects
DATADIR='data/tidy_data'
proj = loadArchRProject(here('data/tidy_data',  'ArchRProjects',
                             'ArchR_DLPFC_multiomeATAC'), showLogo = FALSE)

## grouping cell types by neurogenetic markers
celltypes = getCellColData(proj) %>% as.data.frame() %>%
  arrange(Celltype2) %>% filter(!duplicated(Celltype2)) %>% 
  mutate(cell_class = case_when(grepl('^L[2-6]', Celltype2) ~ 'EXC', 
                                grepl('SST|LAMP|NDN|VIP|PVAL', Celltype2) ~ 'INH',
                                TRUE ~ "GLIA")) %>%
  dplyr::select(Celltype2, cell_class) %>% as.data.frame() %>% 
  split(x = .$Celltype2,f = .$cell_class)
celltypes = unlist(celltypes)

############################################
## get markerPeaks scored and ranked by CNNs
DATADIR='data/tidy_data/celltype_specific_enhancers'
save_top_enh_fn = here(DATADIR, 'rdas', 'rheMac10_DLPFC.candidate_celltype_enhancers.rds')
markerPeaksList = readRDS(save_top_enh_fn) %>% lapply(function(df){ 
  df = df %>% mutate(label =ifelse(is.na(peak2gene),  paste(label, AvgRank, sep = "_"), 
                                   paste(label, AvgRank, peak2gene, sep = "_"))) %>%
    filter(MotifZscore > 0, !is.na(df$peak2gene))
  return(df)
}) %>% lapply(function(df){
  df %>% mutate(
    name = peak,
    peak = gsub('rheMac10:|:250$', '', peak)
  ) %>% dplyr::select(name, peak)%>% deframe() %>%
    GRanges()
})
lengths(markerPeaksList)

FIGDIR='figures/exploratory/celltype_specific_enhancers'
dir.create(here(FIGDIR, 'plots', 'trackplots'), showWarnings = F)
dir.create(here(DATADIR, 'plots'), showWarnings = F)

for(name in names(markerPeaksList)){
  plot_fn = paste0('DLPFC_markerPeaks_SNAIL_candidate_ranked_AllTracks.', name, '.20220217.pdf')
  out_pdf = here(FIGDIR, 'plots',plot_fn)
  # if(!file.exists(out_pdf)){
    plotRegion = GRanges(markerPeaksList[[name]])
    start(plotRegion) = start(plotRegion) - 5e4
    end(plotRegion) = end(plotRegion) + 5e4
    pdf()
    p <- plotBrowserTrack(
      ArchRProj = proj, 
      region = plotRegion,
      groupBy = "Celltype2", 
      useGroups = celltypes,
      title = paste('Marker Peak Tracks for:', name),
      scCellsMax = 2000,
      features = GRangesList(as.list(markerPeaksList[name])),
      loops = getPeak2GeneLinks(proj, corCutOff = 0.45, resolution = 100)
    )
    dev.off()
    plotPDF(p, name = plot_fn, width = 8, height = 8, ArchRProj = proj, addDOC = FALSE)
    system(paste('mv', here('data/tidy_data', 'ArchRProjects','ArchR_DLPFC_multiomeATAC','Plots', plot_fn), 
                 here(FIGDIR, 'plots', 'trackplots')))
  # }
}
