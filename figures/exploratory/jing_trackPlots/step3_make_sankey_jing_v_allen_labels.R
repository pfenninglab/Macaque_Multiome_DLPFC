suppressMessages(library(ArchR))
library(here); library(tidyverse)
library(rtracklayer)
library(future)
library(harmony)

## for plotting 
library(ggsankey)
library(RColorBrewer)

PLOTDIR='figures/exploratory/jing_trackPlots/plots'
in2mm = 25.4

#################################
## 1) grab the scATAC monkey data
proj = loadArchRProject('/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/ArchRProjects/ArchR_DLPFC_scATAC')

## sort the cell types
cellTypes = getCellColData(proj) %>% as.data.frame()%>%
  pivot_longer(cols = c(Celltype2, allen_labels), 
               names_to = 'celltype', values_to = 'values') %>% 
  arrange(desc(grepl('^L[2-3].[A-Z]', values)), desc(grepl('^L[2-3]', values)), 
          desc(grepl('^L[4-5].[A-Z]', values)), desc(grepl('^L[4-5]', values)), 
          desc(grepl('^L6.[A-Z]', values)), desc(grepl('^L6.', values)), 
          desc(grepl('SST|PV', values)),
          desc(grepl('LAMP5|NDNF|VIP|PAX', values)), 
          values) %>%
  pull(values) %>% unique()

cellTypes = cellTypes[c(1:5, 7,6, 8, 10:15, 9, 16:17, 19:20, 18, 21:37 )]

## set the colors for the joint cell types
exc_col = setNames(c(brewer.pal(10, 'Paired'), '#D56DA9'),
                   c('L2.CUX2.MEIS2', 'L4.5.TBX15', 'L4.ALPL', 'L3.CUX2.RORB',
                     'L6.ITGA8','L4.TYR', 'L6.NKD1', 'L6.SYT6', 'L5.PCP4', 
                     'L5.6.NR4A2', 'L5.POU3F1'))
inh_col = setNames(brewer.pal(8, 'Dark2'), cellTypes[21:28])
glia_col = setNames(brewer.pal(9, 'Pastel1'), cellTypes[29:37])

other_cols = cellTypes[! cellTypes %in% (c(exc_col, inh_col, glia_col) %>% names())]
other_cols = setNames(brewer.pal(length(other_cols), 'Set3'), other_cols)
cellTypes_cols = c(exc_col, other_cols, inh_col, glia_col) [cellTypes]
excTypes_cols = cellTypes_cols[1:20]
inhTypes_cols = cellTypes_cols[21:28]

###########################################
## 2) make the sankey diagram of allen to JH clusters, excitatory neurons
df = getCellColData(proj) %>% as.data.frame() %>% as_tibble() 

df_exc_long = df %>% 
  filter(grepl('^L[2-6]', Celltype2)) %>% 
  make_long(Celltype2, allen_labels) %>% 
  mutate(node = factor(node, rev(names(excTypes_cols))))

pdf(here(PLOTDIR, 'macaque_PFC_scATAC_JH_AllenLabels.exc.pdf'), onefile = F, 
    height = 360/in2mm, width = 120/in2mm)
ggplot(df_exc_long, aes(x = x, 
                    next_x = next_x, 
                    node = node, 
                    next_node = next_node,
                    fill = node, 
                    label = node)) +
  geom_sankey() + 
  geom_sankey_label(size = 7, color = "white", fill = "gray40") +
  scale_fill_manual(values = excTypes_cols) +
  labs(x = NULL, y = NULL) + 
  scale_x_discrete(labels=c("Celltype2" = "Rhesus\nMacaque\nLabels", 
                            "allen_labels" = "Allen\nBrain\nLabels"),
                   position = "top") +
  theme_classic(base_size = 26) + 
  theme(legend.position = 'none', 
        axis.line=element_blank(),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())  #remove y axis ticks)
dev.off()


#########################################################################
## 3) make the sankey diagram of allen to JH clusters, excitatory neurons
df_inh_long = df %>% 
  filter( grepl('SST|PV|LAMP5|NDNF|VIP|PAX', Celltype2)) %>% 
  make_long(Celltype2, allen_labels) %>% 
  mutate(node = factor(node, rev(names(inhTypes_cols))))

pdf(here(PLOTDIR, 'macaque_PFC_scATAC_JH_AllenLabels.inh.pdf'), onefile = F, 
    height = 120/in2mm, width = 120/in2mm)
ggplot(df_inh_long, aes(x = x, 
                        next_x = next_x, 
                        node = node, 
                        next_node = next_node,
                        fill = node, 
                        label = node)) +
  geom_sankey() + 
  geom_sankey_label(size = 5, color = "white", fill = "gray40") +
  scale_fill_manual(values = inhTypes_cols) +
  labs(x = NULL, y = NULL) + 
  scale_x_discrete(labels=c("Celltype2" = "Rhesus\nMacaque\nLabels", 
                            "allen_labels" = "Allen\nBrain\nLabels")) +
  theme_void(base_size = 26) + 
  theme(legend.position = 'none', 
        axis.line=element_blank(),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())  #remove y axis ticks)
dev.off()




