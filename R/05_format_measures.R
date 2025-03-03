library(dplyr)
library(tidyr)

scale_fun <- function(v) (v-min(v)) / max(v-min(v))

format_as_dataframe <- function(name, cs){
  # Load the measure results
  load(paste0('measures_', name,'.RData'))
  measures <- as.data.frame(do.call(rbind, measures))

  #we only want the reasonable cluster sizes - set cs
  measures <- measures %>%
    mutate(cluster_size = as.numeric(cluster_size)) %>%
    filter(cluster_size %in% cs)
  
  measures <- measures %>%
    separate(Info, into = c("Name", "Info"), sep = "_&_", extra = "merge", remove = FALSE)
  measures <- measures %>% mutate(Info = sub("_&_", " ", Info))
  measures <- measures %>% mutate(Info = sub('resolution', 'res', Info))
  measures <- measures %>% mutate(Info = sub('cluster_size_', 'cs', Info))
  measures <- measures %>% mutate(
    Name = case_when(
      Name == 'RaceID' ~ 'RID',
      Name == 'PcaKmeans'~ 'PKM',
      Name == 'CIDR'~ 'CIDR',
      Name == 'TsneKmeans'~ 'TKM',
      Name == 'TsneHC'~ 'THC',
      Name == 'PcaHC'~ 'PHC',
      Name == 'pcaReduceS'~ 'PRS',
      Name == 'pcaReduceM'~ 'PRM',
      Name == 'PcaLouvain' ~ 'PLV',
      Name == 'TsneLouvain' ~ 'TLV',
      Name == 'PcaLeiden' ~ 'PLN',
      Name == 'TsneLeiden' ~ 'TLN'
    ))
  
  measures <- measures %>% mutate(
    across(c(DI, IGP, SW, ARI, BSI, BHI), ~ scale_fun(as.numeric(.))), # higher values better
    across(c(AD, ADM, APN, CN), ~ 1 - scale_fun(as.numeric(.))))       # lower values better
  
  save(measures, file = paste0('measures_', name,'_formatted.RData'))
}

format_as_dataframe("BAISE", 2:7)
format_as_dataframe("POLLEN", 2:14)
format_as_dataframe("GOOLAM", 2:8)
format_as_dataframe("DENG", 7:13)
format_as_dataframe("YAN", 5:12)
format_as_dataframe("BEUTNER", 2:30)
