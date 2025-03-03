# REQUIRED LIBRARIES
library(stats)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
library(mclust)
sourceCpp("04_measure_functions/Measures.cpp")

# COMMENTS ARE DESCRIBING BAISE EXAMPLE. THIS IS FOR TESTING.
# md <- list of 2
# md$full <- list of 6
# md$full$cs2 <- classification vector of 49 cells
#
# md$reduced <- list of 49 (one for each cell)
# md$reduced$`cell group 1 removed` <- list of 6
# md$reduced$`cell group 1 removed`$cs2 <- classification vector for 48 cells
# 
# fd <- matrix of cells(columns) and genes(rows)
# 
# index <- list of k-folds to consider for removal
#     e.g. 1:49 would remove 1 cell at a time
#          c(1, 1:48) would remove cell 1 and cell 2 for the first fold,
#           and everything else for the second fold.
#     see the 04_calculate_measures.R and look at the "create_index" function for more insight
#          
# stability = FALSE ; dont do stability measures
# known = FALSE ; dont do external measures
# bd = NULL ; the true classification vector of 49 cells

get_measures <- function(md, fd, index, stability = FALSE, known = FALSE, bd = NULL, verbose=TRUE) {
  
  #initialize measure vector
  measures <- vector();
  
  # calculate the distance matrix of the matrix of cells
  fd.dist <- dist_mat(fd);
  
  # get a vector listing of the nearest neighbor of each cell
  # e.g. cell1 -> 5, 6, 23, 22, 13 means 5 is closest by distance (expression)
  #      followed by 6, etc.
  fd.nn <- list_nn(fd, 50);
  
  #get the name
  method_name <- deparse(substitute(md))
  
  #we now have oc = the list of 6
  oc <- md$full
  
  #we now have rc = the list of 49 cells-removed-lists-of-6
  rc <- md$reduced
  
  #the names of index (cluster sizes)
  params <- names(oc)
  
  #now, we go into the cluster sizes of each. 
  for(p in params) {
    if(verbose) cat('\n>> Working on params ', p,'\n')
    #initialize some holder
    measure_mat <- vector() 
    
    #get the vector list for this cluster size given the full data
    ocp <- oc[[p]]
    
    #rcp is the list of 49 vectors, where each vector is the 2 cluster size of when the i'th cell is removed
    rcp <- lapply(rc, function(x) x[[p]])
    
    #get the cluster size
    cluster_size <- length(unique(ocp))
    
    if(stability){
      #calculate stability measures
      AD <- AverageDistance(ocp = ocp, rcp = rcp, distmat = fd.dist, index = index)
      if(verbose) cat(' AD=',AD)
      ADM <- AverageDistanceBetweenMeans(ocp = ocp, rcp = rcp, fd = fd, index = index)
      if(verbose) cat(' ADM=',ADM)
      APN <- AverageProportionofNonoverlap(ocp = ocp, rcp = rcp, index = index)
      if(verbose) cat(' APN=',APN)
      #record stability measures
      measure_mat <- c(measure_mat, AD=AD, ADM=ADM, APN=APN)
    }
    
    if(known){
      ARI <- mclust::adjustedRandIndex(ocp, bd)
      if(verbose) cat(' ARI=',ARI)
      measure_mat <- c(measure_mat, ARI=ARI)
      
      if(stability){
        BHI <- BiologicalHomogeneityIndex(ocp, bd) 
        if(verbose) cat(' BHI=',BHI)
        BSI <- BiologicalStabilityIndex(ocp, bd, rcp )
        if(verbose) cat(' BSI=',BSI)
        measure_mat <- c(measure_mat, BHI=BHI, BSI=BSI)
      }
    }
    
    CN <- Connectivity(ocp = ocp, nn = fd.nn, h = 5)
    if(verbose) cat(' CN=',CN)
    DI <- DunnIndexCpp(ocp = ocp, distmat = fd.dist)
    if(verbose) cat(' DI=',DI)
    IGP <- InGroupProportion(ocp = ocp, nn = fd.nn)
    if(verbose) cat(' IGP=',IGP)
    SW <- SilhouetteDistancecpp(ocp = ocp, distmat = fd.dist)
    if(verbose) cat(' SW=',SW)
    if(is.na(SW)) message('\n Only 1 cluster. SW cannot be calculated. Poor choice in parameters.')
    measure_mat <- c(measure_mat, CN=CN, DI=DI, IGP=IGP, SW=SW)
    #update full list of measures
    measure_mat <- c(cluster_size = cluster_size, measure_mat)
    #compile this round of measures
    measures <- rbind(measures, measure_mat);
  }
  rownames(measures) <- params
  return(measures)
}

AverageDistance <- function (ocp, rcp, distmat, index) {
  score = 0; 
  n = length(index)
  
  for ( i in 1:n ) {  
    rcjk <- rcp[[ i ]];
    ocjk <- ocp[ -index[[ i ]] ]
    njk <- length(ocjk) 
    tab <- as.data.frame(table(rcjk, ocjk)) 
    tab <- tab[tab[,3] != 0,] 
    
    new <- t(apply( tab , 1 , function(x) { 
        x <- as.numeric(x) 
        rcjk_ind <- which(rcjk == x[1])
        ocjk_ind <-  which(ocjk == x[2]) 
        value <- x[3] * mean(distmat[rcjk_ind,ocjk_ind])
        
        return(value) 
        }))
    score = score + sum(new) / njk
  }
  score = score / n 
  return(score)
}


AverageDistanceBetweenMeans <- function (ocp, rcp, fd, index) {
  score = 0;
  n = length(index)
  fd_ind <- 1:length(ocp)
  
  for ( i in 1:n ) {
    rcjk <- rcp[[ i ]];
    ocjk <- ocp[ -index[[ i ]] ];
    njk <- length(ocjk)
    
    fdjk_ind <- fd_ind[ -index[[i]] ]
    
    tab <- as.data.frame(table(rcjk, ocjk))
    tab <- tab[tab[,3] != 0,]
    tab <- t(apply( tab , 1 , function(x) {
      x <- as.numeric(x)
      
      rcjk_ind <- fdjk_ind[which(rcjk == x[1])]
      ocjk_ind <- fdjk_ind[which(ocjk == x[2])]
      rcjk_cent <- rowMeansC(as.matrix(fd[,rcjk_ind]))
      ocjk_cent <- rowMeansC(as.matrix(fd[,ocjk_ind]))
      x[3] <- x[3] * euclideanDistance(rcjk_cent, ocjk_cent)
      
      return(x)
      }))
    score = score + sum(tab[,3]) / njk
  }
  score = score / n
  return(score)
}


AverageProportionofNonoverlap <- function (ocp, rcp, index) {
  score = 0;
  n = length(index)
  
  for ( i in 1:n ) {
    rcjk <- rcp[[ i ]];
    ocjk <- ocp[ -index[[ i ]] ];
    njk <- length(ocjk)
    tab <- as.data.frame(table(rcjk, ocjk))
    tab <- tab[tab[,3] != 0,]
    tab <- t(apply( tab , 1 , function(x) {
      x <- as.numeric(x)
      
      numer <- length( intersect( which(rcjk == x[1]), which(ocjk == x[2]) ) )
      denom <- length( which(ocjk == x[2]) )
      x[3] <- x[3] * (1 - numer/denom)
      
      return(x)
    }))
    
    score = score + sum(tab[,3]) / njk
  }
  score = score / n
  return(score)
}




BiologicalHomogeneityIndex <- function(oc, bd){
  uck = unique(oc);
  k = length(uck);
  res = 0;
  
  for(i in 1:k){
    oc_index = which(oc == uck[i]);
    ni <- length(oc_index);
    m=0;
    if(ni > 1){
      for(j in 1:(ni-1)) for(h in (j+1):ni){ # looping j < h, which is 0.5 * j!=h
        if(bd[oc_index[j]] == bd[oc_index[h]]) m = m + 1;
      }
      res = res + m / (ni*(ni-1)/2);  
    }
  } 
  return(res / k)
}

BiologicalStabilityIndex <- function(oc, bd, rc){
  bd <- as.numeric(as.factor(bd))
  sol <- BioSIndexCpp(oc, bd, rc)
  return(sol);
}


Connectivity <- function(ocp, nn, h) {
  score <- 0
  
  for ( i in 1:length(ocp) ) {
    check <- ocp[i] != ocp[nn[[i]][1:h]]
    if (any(check)) {
      score = score + sum(1/which(check))
    }
  }
  
  return(score)
}

# DunnIndex <- function(ocp, distmat) {
#   ucp <- unique(ocp)
#   m <- length(ucp); if(m == 1) return(0);
#   num_dist <- c();
#   for(i in 1:(m-1)){
#     for(j in (i+1):m){
#       sc_index = which(ocp == ucp[i]);
#       tc_index = which(ocp == ucp[j]);
#       num_dist <- c(num_dist, min(distmat[sc_index, tc_index]));
#     }
#   }
#   num <- min(num_dist);
#   denom_dist <- c();
#   for(i in 1:m){
#     rc_index = which(ocp == ucp[i]);
#     denom_dist <- c(denom_dist, max(distmat[rc_index, rc_index]))
#   }
#   denom <- max(denom_dist);
#   score <- num / denom;
#   return(score);
# }


InGroupProportion <- function(ocp, nn) {
  ucp <- unique(ocp)
  m <- length(ucp)
  
  score = 0;
  for ( i in 1:m ) {
    ocp_index = which(ocp == ucp[i]);
    score <- score + (sum(unlist(lapply(ocp_index, function(x){
      nn[[x]][1]
      })) %in% ocp_index) / length(ocp_index));
  }
  return(score/m);
}






