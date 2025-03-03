# install packages
#
#install.packages("BiocManager");
library(BiocManager)
#BiocManager::install('pcaMethods', update = FALSE)
library(pcaMethods)
#install.packages('mnormt')
library(mnormt)
#install.packages("mclust")
library(mclust)
#install.packages("clue")
library(clue)

#setwd('/home/owenvisser/ClustPaper')
#install.packages("pcaReduce", repos = NULL, type="source")
library(pcaReduce)


# pca reduce method
runpcaReduce <- function(data, cs, ...){
  top <- max(max(cs) - 1, 2) ; bot <- max(min(cs),2)
  
  if(top == 2 & bot == 2){ #niche case hard-coded
    index <- 2; sizeIndex <- (top+1):bot
  }else{
    index <- length(cs):1; sizeIndex <- (top+1):bot
  }

  bycluster = list();
  listOfResults = PCAreduce(D_t = t(data), q = top, ...)
  for(i in index){ 
    message('working on the following cluster size: ', sizeIndex[i])
    parts <- lapply(listOfResults, function(x) as.cl_partition(x[,i])) 
    combined <- cl_consensus(parts, method = "HE", control = list(nruns = 50, k = sizeIndex[i]))
    assignments <- apply(combined$.Data, 1, function(row) which(row == max(row)))
    bycluster[[paste0('cs', sizeIndex[i])]] <- assignments
  }
  return(bycluster)
}

# 
# setwd('/home/owenvisser/ClustPaper')
# 
# # import data
# load('m3m.RData')
# 
# ################################################## method = 'M'
# full <-  runpcaReduce(data = info_list$data, cs = info_list$cluster_sizes, nbt = 100, method = 'M')
# 
# reduced <- list();
# cntr <- 0
# for (i in info_list$index) {
#   cntr = cntr + 1
#   message('running reduced set ', cntr)
#   
#   hold <- runpcaReduce(d = info_list$data[,-i], cs = info_list$cluster_sizes, nbt = 100, method = 'M')
#   reduced[[paste('cell group ', cntr, ' removed')]] <- hold
# }
# 
# pcaReduceM <- list(reduced = reduced, full = full)
# save(pcaReduceM,  file = 'm3m_pcaReduceM.RData')
# 
# ################################################## method = 'S'
# full <-  runpcaReduce(data = info_list$data, cs = info_list$cluster_sizes, nbt = 100, method = 'S')
# 
# reduced <- list();
# cntr <- 0
# for (i in info_list$index) {
#   cntr = cntr + 1
#   message('running reduced set ', cntr)
#   
#   hold <- runpcaReduce(d = info_list$data[,-i], cs = info_list$cluster_sizes, nbt = 100, method = 'S')
#   reduced[[paste('cell group ', cntr, ' removed')]] <- hold
# }
# 
# pcaReduceS <- list(reduced = reduced, full = full)
# save(pcaReduceS,  file = 'm3m_pcaReduceS.RData')
# 
# # clear cache
# rm(list = ls())

