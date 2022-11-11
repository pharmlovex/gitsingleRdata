
## Cell level filtering 
merge_filtered <- subset(merge_obj, subset = nUMI > 1000 &
                           nUMI < 15000 &
                           nGene > 500 &
                           mitoPercent > 0.5 & 
                           log10GenesPerUMI > 0.85 &
                           riboPercent >1 &
                           riboPercent < 25)