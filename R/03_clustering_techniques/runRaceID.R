# install packages
# 
#install.packages('RaceID')
library(RaceID)
library(Matrix)

# raceID function
# 
runRaceID <- function(d, cs){
  scm <- SCseq(d)
  
  #null filtering, no pre-processing will be done.
  scm <- filterdata(scm, mintotal = 1, minexpr = 1,
                    minnumber = 1, bmode = "RaceID", verbose = FALSE)
  
  scm <- compdist(scm, metric = "pearson", FSelect = FALSE)
  
  bycluster = list();
  for(i in cs){
    message('working on ', i)
    hold <- clustexp(scm, sat = FALSE, cln = i,
                     clustnr = 30, bootnr = 50, rseed = 17000,
                     FUNcluster = "kmedoids", verbose = FALSE)
    bycluster[[paste0('cs',i)]] <- hold@cluster$kpart
  }
  
  return(bycluster)
}

# setwd('/home/owenvisser/ClustPaper')
# 
# # import data
# load('m3m.RData')
# 
# # get full cluster
# full <-  runRaceID(d = info_list$data,
#                  cs = info_list$cluster_sizes)
# 
# # get reduced clusters
# reduced <- list();
# cntr <- 0
# for (i in info_list$index) {
#   cntr = cntr + 1
#   message('running reduced set ', cntr)
#   hold <- runRaceID(d = info_list$data[, -i],
#                   cs = info_list$cluster_sizes)
#   reduced[[paste('cell group ', cntr, ' removed')]] <- hold
# }
# 
# # compile results
# RaceID <- list(reduced = reduced, full = full)
# 
# # save results
# save(RaceID,  file = 'm3m_RaceID.RData')
# 
# # clear cache
# rm(list = ls())