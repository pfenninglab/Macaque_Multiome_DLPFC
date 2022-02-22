
## get 
celltypes1 = c("EXC","INH_LAMP5","INH_PVALB","INH_SST" ,
               "INH_VIP", "Astro", "Endo", "Microglia","Oligo", "OPC")
celltypes2 = c("L2.CUX2.MEIS2", "L3.CUX2.RORB" ,  "L4.5.TBX15" , "L4.ALPL" ,
               "L4.TYR" , "L5.6.NR4A2" ,"L5.PCP4" ,  "L5.POU3F1", "L6.ITGA8", 
               "L6.NKD1",  "L6.SYT6", 
               "LAMP5",  "NDNF",  "VIP",  "PV.BC" , "PV.ChC" ,"SST", 
               "Astro", "Endo" , "Microglia", "Mural" , "Oligo" , "OPC" )
celltypes2_cols = setNames( c(brewer.pal(11, 'Paired'), brewer.pal(6, 'Pastel1'), 
                              brewer.pal(6, 'Dark2')), celltypes2)

cellsList = lapply(model_list, function(model){
  tmp = apply(model$h, 1, function(x) x / max(x, na.rm = T)) %>% as.data.frame()
  mycols = names(tmp) = gsub('V', 'K', names(tmp))
  df = colData(peakMat) %>% as.data.frame() %>%
    dplyr::select(Sample, Region, Celltype1, Celltype2)
  df = cbind(df, tmp) %>%
    pivot_longer(cols = starts_with('K'), values_to = 'loadings', names_to = 'factor') %>%
    group_by(Celltype2, factor) %>% summarise(loadings = mean(loadings)) %>%
    ungroup() %>%  
    mutate(factor = factor(factor, levels = names(tmp)), 
           Celltype2 = factor(Celltype2, celltypes2))
})

pdf('tmp.pdf', height = 8, width = 12)
ggplot(df, aes(x = factor, y = loadings, fill = Celltype2)) + 
  geom_bar(stat = 'identity', position="fill") + 
  scale_fill_manual(values = celltypes2_cols) + 
  theme_bw() 
dev.off()


