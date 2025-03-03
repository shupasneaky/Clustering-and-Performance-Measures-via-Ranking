# Load all techniques - may need to manually set directory
for(f in list.files("03_clustering_techniques")) source(paste0("03_clustering_techniques/",f))

## for CIDR
do_CIDR <- function(name){
  full = runCIDR(d=sce@assays@data$logcounts, cs = cluster_sizes)
  reduced = list()
  for(i in 1:ncol(sce)){
    message('Processing cell ',i)
    reduced[[paste0('cell ', i, 'removed')]] <- runCIDR(d=sce@assays@data$logcounts[,-i], cs = cluster_sizes)
  }
  CIDR <- list(reduced = reduced, full = full)
  save(CIDR,  file = paste0(name,'_CIDR.RData'))
  message(name,' processing complete for CIDR\n','File saved as ',paste0(name,'_CIDR.RData'))
  gc()
}

## for pcaReduce
do_pcaReduce <- function(name){
  # M method
  full <- runpcaReduce(data = sce@assays@data$logcounts, cs = cluster_sizes, nbt = 100, method = 'M')
  reduced = list()
  for(i in 1:ncol(sce)){
    message('Processing cell ',i)
    reduced[[paste0('cell ', i, 'removed')]] <- runpcaReduce(data = sce@assays@data$logcounts[,-i], cs = cluster_sizes, nbt = 100, method = 'M')
  }
  pcaReduceM <- list(reduced = reduced, full = full)
  save(pcaReduceM,  file = paste0(name,'_pcaReduceM.RData'))
  message(name,' processing complete for pcaReduceM\n','File saved as ',paste0(name,'_pcaReduceM.RData'))
  gc()
  
  # S method
  full <- runpcaReduce(data = sce@assays@data$logcounts, cs = cluster_sizes, nbt = 100, method = 'S')
  reduced = list()
  for(i in 1:ncol(sce)){
    message('Processing cell ',i)
    reduced[[paste0('cell ', i, 'removed')]] <- runpcaReduce(data = sce@assays@data$logcounts[,-i], cs = cluster_sizes, nbt = 100, method = 'S')
  }
  pcaReduceS <- list(reduced = reduced, full = full)
  save(pcaReduceS,  file = paste0(name,'_pcaReduceS.RData'))
  message(name,' processing complete for pcaReduceS\n','File saved as ',paste0(name,'_pcaReduceS.RData'))
  gc()
}

do_RaceID <- function(name){
  full <- runRaceID(d = sce@assays@data$logcounts, cs = cluster_sizes)
  reduced = list()
  for(i in 1:ncol(sce)){
    message('Processing cell ',i)
    reduced[[paste0('cell ', i, 'removed')]] <- runRaceID(d = sce@assays@data$logcounts[,-i], cs = cluster_sizes)
  }
  RaceID <- list(reduced = reduced, full = full)
  save(RaceID,  file = paste0(name,'_RaceID.RData'))
  message(name,' processing complete for RaceID\n','File saved as ',paste0(name,'_RaceID.RData'))
  gc()
}

do_Homebrew <- function(name){
  d <- get.dr(d = sce@assays@data$logcounts)
  full.pkm <- runPcaKmeans(pca = d$pca, cs = cluster_sizes)
  full.tkm <- runTsneKmeans(tsne = d$tsne, cs = cluster_sizes)
  full.phc <- runPcaHC(dmat = d$pca.d, cs = cluster_sizes)
  full.thc <- runTsneHC(dmat = d$tsne.d, cs = cluster_sizes)
  full.plv <- runPcaLouvain(dmat = d$pca.dmat, resolution = resolution, knn = knn)
  full.tlv <- runTsneLouvain(dmat = d$tsne.dmat, resolution = resolution, knn = knn)
  full.pln <- runPcaLeiden(dmat = d$pca.dmat, resolution = resolution, knn = knn)
  full.tln <- runTsneLeiden(dmat = d$tsne.dmat, resolution = resolution, knn = knn)
  
  # get reduced clusters
  reduced.pkm <- list();
  reduced.tkm <- list();
  reduced.phc <- list();
  reduced.thc <- list();
  reduced.plv <- list();
  reduced.tlv <- list();
  reduced.pln <- list();
  reduced.tln <- list();
  
  for (i in 1:ncol(sce)) {
    message('Processing cell ',i, ' for all homebrew techniques...')
    d <- get.dr(d = sce@assays@data$logcounts[,-i])
    reduced.pkm[[paste('cell group ', i, ' removed')]] <- runPcaKmeans(pca = d$pca, cs = cluster_sizes)
    reduced.tkm[[paste('cell group ', i, ' removed')]] <- runTsneKmeans(tsne = d$tsne, cs = cluster_sizes)
    reduced.phc[[paste('cell group ', i, ' removed')]] <- runPcaHC(dmat = d$pca.d, cs = cluster_sizes)
    reduced.thc[[paste('cell group ', i, ' removed')]] <- runTsneHC(dmat = d$tsne.d, cs = cluster_sizes)
    reduced.plv[[paste('cell group ', i, ' removed')]] <- runPcaLouvain(dmat = d$pca.dmat, resolution = resolution, knn = knn)
    reduced.tlv[[paste('cell group ', i, ' removed')]] <- runTsneLouvain(dmat = d$tsne.dmat, resolution = resolution, knn = knn)
    reduced.pln[[paste('cell group ', i, ' removed')]] <- runPcaLeiden(dmat = d$pca.dmat, resolution = resolution, knn = knn)
    reduced.tln[[paste('cell group ', i, ' removed')]] <- runTsneLeiden(dmat = d$tsne.dmat, resolution = resolution, knn = knn)
  }
  
  # compile results
  PcaKmeans <- list(reduced = reduced.pkm, full = full.pkm)
  TsneKmeans <- list(reduced = reduced.tkm, full = full.tkm)
  PcaHC <- list(reduced = reduced.phc, full = full.phc)
  TsneHC <- list(reduced = reduced.thc, full = full.thc)
  PcaLouvain <- list(reduced = reduced.plv, full = full.plv)
  TsneLouvain <- list(reduced = reduced.tlv, full = full.tlv)
  PcaLeiden <- list(reduced = reduced.pln, full = full.pln)
  TsneLeiden <- list(reduced = reduced.tln, full = full.tln)
  
  # save results
  save(PcaKmeans,  file = paste0(name,'_PcaKmeans.RData'))
  save(TsneKmeans,  file = paste0(name,'_TsneKmeans.RData'))
  save(PcaHC,  file = paste0(name,'_PcaHC.RData'))
  save(TsneHC,  file = paste0(name,'_TsneHC.RData'))
  save(PcaLouvain,  file = paste0(name,'_PcaLouvain.RData'))
  save(TsneLouvain,  file = paste0(name,'_TsneLouvain.RData'))
  save(PcaLeiden,  file = paste0(name,'_PcaLeiden.RData'))
  save(TsneLeiden,  file = paste0(name,'_TsneLeiden.RData'))
  gc()
}

#load("BAISE_data_preprocessed.RData")
load("Baise_sc_data/BAISE_data_preprocessed.RData")
cluster_sizes = 2:7
resolution = 1:9 / 10
knn = c(10, 25, 40)
do_CIDR("BAISE")
do_pcaReduce("BAISE")
do_RaceID("BAISE")
do_Homebrew("BAISE")

load("Deng_sc_data/DENG_data_preprocessed.RData")
cluster_sizes = 7:13
resolution = 1:9 / 10
knn = c(10, 25, 40)
do_CIDR("DENG")
do_pcaReduce("DENG")
do_RaceID("DENG")
do_Homebrew("DENG")

load("Goolam_sc_data/GOOLAM_data_preprocessed.RData")
cluster_sizes = 2:8
resolution = 1:9 / 10
knn = c(10, 25, 40)
do_CIDR("GOOLAM")
do_pcaReduce("GOOLAM")
do_RaceID("GOOLAM")
do_Homebrew("GOOLAM")

load("Yan_sc_data/YAN_data_preprocessed.RData")
cluster_sizes = 5:12
resolution = 1:9 / 10
knn = c(10, 25, 40)
do_CIDR("YAN")
do_pcaReduce("YAN")
do_RaceID("YAN")
do_Homebrew("YAN")

load("Pollen_sc_data/POLLEN_data_preprocessed.RData")
cluster_sizes = 2:14
resolution = 1:9 / 10
knn = c(10, 25, 40)
do_CIDR("POLLEN")
do_pcaReduce("POLLEN")
do_RaceID("POLLEN")
do_Homebrew("POLLEN")

load("Beutner_mice_age_sc_data/BEUTNER_data_preprocessed.RData")
cluster_sizes = 2:30
resolution = 1:9 / 10
knn = c(20, 30, 40, 50)
do_CIDR("BEUTNER")
do_pcaReduce("BEUTNER")
do_RaceID("BEUTNER")
do_Homebrew("BEUTNER")
