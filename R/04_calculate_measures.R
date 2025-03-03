# Load in all the results...
# setwd("C:/Users/owvis/UFL Dropbox/Owen Visser/Datasets/Beutner_mice_age_sc_data")
# setwd('/home/owenvisser/ClustPaper')

source('04_measure_functions/Measures.R')


do_measures <- function(name, stab=TRUE, knwn=TRUE, fold_index = NULL, knwn_vec = NULL){
  files = list.files(pattern = name)
  for(f in files) load(f)
  
  if(is.null(fold_index)) stop('ensure fold index is described')
  if(is.null(knwn_vec) && knwn==T) stop('can\'t calculate known if you dont supply the true vector')
  
  # Get the list of methods so we don't have to do anything manually
  # Can always just define it manually if the list of methods is not very long
  methods <- files[grepl(paste0(name,'_'), files)]
  methods <- gsub(paste0(name,'_'), '', methods)
  methods <- gsub('.RData', '', methods)
  
  # Below will now reference the first method in 'methods', and you can subset it
  # get(methods[1])$full
  # source("C:/Users/owvis/UFL Dropbox/Owen Visser/Dissertation/ClustPaper/R/Measures_Set.R")
  
  measures <- list();
  for (m in methods) {
    cat('\n\n ===Method ', m,'===')
    measures[[m]] <- get_measures(md = get(m),
                                  fd = sce@assays@data$logcounts,
                                  index = fold_index,
                                  stability = stab,
                                  known = knwn,
                                  bd = knwn_vec)
    measures[[m]] <- cbind(Info = paste0(m,'_&_', rownames(measures[[m]])), measures[[m]])
    rownames(measures[[m]]) <- NULL
  }
  save(measures, file = paste0('measures_',name,'.RData'))
  message('\n\n Methods complete. File saved as ',paste0('measures_',name,'.RData'))
}

create_index <- function(len, k) {
  sample_vec <- sample.int(len)
  sample_sizes <- floor(len/k)
  sr <- split(sample_vec, cut(seq_along(1:len), breaks = k, labels = FALSE))
  return(sr)
}

load("Baise_sc_data/BAISE_data_preprocessed.RData")
do_measures("BAISE", fold_index = create_index(len=ncol(sce), k=ncol(sce)), knwn_vec = gsub("_.*", "", colnames(sce)))

load("Deng_sc_data/DENG_data_preprocessed.RData")
do_measures("DENG", fold_index = create_index(len=ncol(sce), k=ncol(sce)), knwn_vec = gsub("_.*", "", colnames(sce)))

load("Goolam_sc_data/GOOLAM_data_preprocessed.RData")
do_measures("GOOLAM", fold_index = create_index(len=ncol(sce), k=ncol(sce)), knwn_vec = gsub("_.*", "", colnames(sce)))

load("Yan_sc_data/YAN_data_preprocessed.RData")
do_measures("YAN", fold_index = create_index(len=ncol(sce), k=ncol(sce)), knwn_vec = gsub("_.*", "", colnames(sce)))

load("Pollen_sc_data/POLLEN_data_preprocessed.RData")
do_measures("POLLEN", fold_index = create_index(len=ncol(sce), k=ncol(sce)), knwn_vec = gsub("_.*", "", colnames(sce)))

load("Beutner_mice_age_sc_data/BEUTNER_data_preprocessed.RData")
do_measures("BEUTNER", fold_index = create_index(len=ncol(sce), k=100), knwn_vec = gsub("_.*", "", colnames(sce)))

