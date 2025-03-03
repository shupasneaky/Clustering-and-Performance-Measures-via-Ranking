#Run Full

#libraries
library(cccd)
library(Rtsne)
library(stats)
library(igraph)

get.dr <- function(d){
  pca <- prcomp(x = t(d), rank. = 50)$x
  pca.d <- dist(pca, method = "euclidean")
  pca.dmat <- as.matrix(pca.d)
  
  tsne <- Rtsne(
    pca.dmat,
    dims = 3,
    perplexity = (ncol(pca.dmat)-1)/3,
    theta = 0.5,
    check_duplicates = FALSE,
    pca = FALSE,
    is_distance = TRUE,
    Y_init = NULL,
    normalize = FALSE
  )$Y
  tsne.d <- dist(tsne, method = "euclidean")
  tsne.dmat <- as.matrix(tsne.d)
  
  dimnames(pca) <- dimnames(pca.d) <- dimnames(pca.dmat) <- NULL
  dimnames(tsne) <- dimnames(tsne.d) <- dimnames(tsne.dmat) <- NULL
  
  return(list(pca = pca, pca.d = pca.d, pca.dmat = pca.dmat,
              tsne = tsne, tsne.d = tsne.d, tsne.dmat = tsne.dmat))
}

runPcaKmeans <- function(pca, cs){
  
  bycluster = list();
  for (i in cs) {
    message('working on ', i)
    km <- kmeans(pca, centers = i, iter.max = 100, nstart = 10, trace=FALSE)
    bycluster[[paste0('cs',i)]] <- km$cluster
  }
  return(bycluster)
}

runTsneKmeans <- function(tsne, cs){
  bycluster = list();
  for (i in cs) {
    message('working on ', i)
    km <- kmeans(tsne, centers = i, iter.max = 100, nstart = 10, trace=FALSE)
    bycluster[[paste0('cs',i)]] <- km$cluster
  }
  return(bycluster)
}

runPcaHC <- function(dmat, cs){
  bycluster = list();
  for (i in cs) {
    message('working on ', i)
    tree <- hclust(d = dmat, method = "ward.D")
    bycluster[[paste0('cs',i)]] <- cutree(tree = tree, k = i)
  }
  return(bycluster)
}

runTsneHC <- function(dmat, cs){
  bycluster = list();
  for (i in cs) {
    message('working on ', i)
    tree <- hclust(d = dmat, method = "ward.D")
    bycluster[[paste0('cs',i)]] <- cutree(tree = tree, k = i)
  }
  return(bycluster)
}

runPcaLouvain <- function(dmat, resolution, knn){
  bycluster = list();
  for(j in knn) {
    graph <- nng(dx = dmat, mutual = TRUE, k = j)
    for (i in resolution) {
      message('working on ', j, ' ', i)
      h <- cluster_louvain(graph, resolution = i)
      bycluster[[paste0('knn_',j,'_&_resolution_',i)]] <- h$membership
    }
  }
  return(bycluster)
}

runTsneLouvain <- function(dmat, resolution, knn){
  bycluster = list();
  for(j in knn) {
    graph <- nng(dx = dmat, mutual = TRUE, k = j)
    for (i in resolution) {
      message('working on ', j, ' ', i)
      h <- cluster_louvain(graph, resolution = i)
      bycluster[[paste0('knn_',j,'_&_resolution_',i)]] <- h$membership
    }
  }
  return(bycluster)
}


runPcaLeiden <- function(dmat, resolution, knn){
  bycluster = list();
  for(j in knn) {
    graph <- nng(dx = dmat, mutual = TRUE, k = j)
    for (i in resolution) {
      message('working on ', j, ' ', i)
      h <- cluster_leiden(graph, resolution_parameter = i)
      bycluster[[paste0('knn_',j,'_&_resolution_',i)]] <- h$membership
    }
  }
  return(bycluster)
}

runTsneLeiden <- function(dmat, resolution, knn){
  bycluster = list();
  for(j in knn) {
    graph <- nng(dx = dmat, mutual = TRUE, k = j)
    for (i in resolution) {
      message('working on ', j, ' ', i)
      h <- cluster_leiden(graph, resolution_parameter = i)
      bycluster[[paste0('knn_',j,'_&_resolution_',i)]] <- h$membership
    }
  }
  return(bycluster)
}
# 
# setwd('/home/owenvisser/ClustPaper')
# 
# # import data
# load('m3m.RData')
# 
# d <- get.dr(d = info_list$data)
# full.pkm <- runPcaKmeans(pca = d$pca, cs = info_list$cluster_sizes)
# full.tkm <- runTsneKmeans(tsne = d$tsne, cs = info_list$cluster_sizes)
# full.phc <- runPcaHC(dmat = d$pca.d, cs = info_list$cluster_sizes)
# full.thc <- runTsneHC(dmat = d$tsne.d, cs = info_list$cluster_sizes)
# full.plv <- runPcaLouvain(dmat = d$pca.dmat, resolution = info_list$resolution, knn = info_list$knn)
# full.tlv <- runTsneLouvain(dmat = d$tsne.dmat, resolution = info_list$resolution, knn = info_list$knn)
# full.pln <- runPcaLeiden(dmat = d$pca.dmat, resolution = info_list$resolution, knn = info_list$knn)
# full.tln <- runTsneLeiden(dmat = d$tsne.dmat, resolution = info_list$resolution, knn = info_list$knn)
# 
# # get reduced clusters 
# reduced.pkm <- list();
# reduced.tkm <- list();
# reduced.phc <- list();
# reduced.thc <- list();
# reduced.plv <- list();
# reduced.tlv <- list();
# reduced.pln <- list();
# reduced.tln <- list();
# 
# cntr <- 0
# for (i in info_list$index) {
#   cntr = cntr + 1
#   message('running reduced set ', cntr)
#   d <- get.dr(d = info_list$data[, -i])
#   reduced.pkm[[paste('cell group ', cntr, ' removed')]] <- runPcaKmeans(pca = d$pca, cs = info_list$cluster_sizes)
#   reduced.tkm[[paste('cell group ', cntr, ' removed')]] <- runTsneKmeans(tsne = d$tsne, cs = info_list$cluster_sizes)
#   reduced.phc[[paste('cell group ', cntr, ' removed')]] <- runPcaHC(dmat = d$pca.d, cs = info_list$cluster_sizes)
#   reduced.thc[[paste('cell group ', cntr, ' removed')]] <- runTsneHC(dmat = d$tsne.d, cs = info_list$cluster_sizes)
#   reduced.plv[[paste('cell group ', cntr, ' removed')]] <- runPcaLouvain(dmat = d$pca.dmat, resolution = info_list$resolution, knn = info_list$knn)
#   reduced.tlv[[paste('cell group ', cntr, ' removed')]] <- runTsneLouvain(dmat = d$tsne.dmat, resolution = info_list$resolution, knn = info_list$knn)
#   reduced.pln[[paste('cell group ', cntr, ' removed')]] <- runPcaLeiden(dmat = d$pca.dmat, resolution = info_list$resolution, knn = info_list$knn)
#   reduced.tln[[paste('cell group ', cntr, ' removed')]] <- runTsneLeiden(dmat = d$tsne.dmat, resolution = info_list$resolution, knn = info_list$knn)
# }
# 
# # compile results
# PcaKmeans <- list(reduced = reduced.pkm, full = full.pkm)
# TsneKmeans <- list(reduced = reduced.tkm, full = full.tkm)
# PcaHC <- list(reduced = reduced.phc, full = full.phc)
# TsneHC <- list(reduced = reduced.thc, full = full.thc)
# PcaLouvain <- list(reduced = reduced.plv, full = full.plv)
# TsneLouvain <- list(reduced = reduced.tlv, full = full.tlv)
# PcaLeiden <- list(reduced = reduced.pln, full = full.pln)
# TsneLeiden <- list(reduced = reduced.tln, full = full.tln)
# 
# # save results
# save(PcaKmeans,  file = 'm3m_PcaKmeans.RData')
# save(TsneKmeans,  file = 'm3m_TsneKmeans.RData')
# save(PcaHC,  file = 'm3m_PcaHC.RData')
# save(TsneHC,  file = 'm3m_TsneHC.RData')
# save(PcaLouvain,  file = 'm3m_PcaLouvain.RData')
# save(TsneLouvain,  file = 'm3m_TsneLouvain.RData')
# save(PcaLeiden,  file = 'm3m_PcaLeiden.RData')
# save(TsneLeiden,  file = 'm3m_TsneLeiden.RData')
# 
# # clear cache
# rm(list = ls())
