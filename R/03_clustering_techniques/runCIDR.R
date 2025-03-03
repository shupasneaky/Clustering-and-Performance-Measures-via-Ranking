## CIDR
######################################
library(devtools)
#devtools::install_github("VCCRI/CIDR", upgrade = 'never')
library(cidr)
library(stats)


runCIDR <- function(d, cs){
  scd <- scDataConstructor(d)
  scd <- determineDropoutCandidates(scd)
  scd <- wThreshold(scd)
  scd <- scDissim(scd, threads = 1)
  scd <- scPCA(scd, plotPC = FALSE)
  scd <- nPC(scd)
  
  bycluster = list();
  for(i in cs){
    message('working on ', i)
    hold <- scCluster(scd, nCluster = i, nPC = scd@nPC, cMethod = "ward.D2")
    bycluster[[paste0('cs',i)]] <- hold@clusters
  }
  return(bycluster)
}

# # 
# setwd('/home/owenvisser/ClustPaper')
# 
# # import data
# load('m3m.RData')
# 
# # get full cluster
# full <-  runCIDR(d = info_list$data,
#                  cs = info_list$cluster_sizes)
# 
# # get reduced clusters
# reduced <- list();
# cntr <- 0
# for (i in info_list$index) {
#   cntr = cntr + 1
#   message('running reduced set ', cntr)
#   hold <- runCIDR(d = info_list$data[, -i],
#                   cs = info_list$cluster_sizes)
#   reduced[[paste('cell group ', cntr, ' removed')]] <- hold
# }
# 
# # compile results
# CIDR <- list(reduced = reduced, full = full)
# 
# # save results
# save(CIDR,  file = 'm3m_CIDR.RData')
# 
# # clear cache
# rm(list = ls())
# 