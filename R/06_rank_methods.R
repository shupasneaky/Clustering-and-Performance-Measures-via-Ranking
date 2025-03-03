##Packages
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(xtable)
library(RankAggreg)

## Run rank aggregation on the measures. 
getRA <- function(d,
                  imp = NULL,
                  m = NULL){
  
  if (is.null(m)) stop('must select measure names')
  if (is.null(imp)) imp=rep(1,length(m))
  
  sizes <- unique(d$cluster_size)
  size_list <- list();
  
  for(s in sizes){
    h <- d %>% filter(cluster_size == s) %>%
      unite('all', Name:cluster_size, remove = TRUE, sep='_&_' ) %>%
      column_to_rownames(var = 'all') %>% t()
    
    if(ncol(h) < 2){
      size_list[[s]] <- colnames(h)
    }else{
      h_list <- t(apply(h, 1, function(x){ colnames(h)[order(x, decreasing = TRUE)] }))
      h_weights <- t(apply(h, 1, function(x){ x[order(x, decreasing = TRUE)] }))
      h_ranks <- RankAggreg(h_list,
                            seed = 0,
                            weights = h_weights,
                            k = min(ncol(h_list), 10), # have at least 10 in the ranking to ensure consistency.
                            importance = imp, # equal to default unless manually set.
                            standardizeWeights = FALSE, # already standardized
                            verbose=FALSE) # no talky
      
      size_list[[paste(s)]] <- h_ranks$top.list[1:min(ncol(h_list), 1)]
    }
  }
  
  top_d <- d %>% unite('all', Name:cluster_size, remove = TRUE, sep='_&_' ) %>%
    filter(all %in% unlist(size_list)) %>%
    column_to_rownames(var = 'all') %>%
    t()
  d_lists <- t(apply(top_d, 1, function(x){ colnames(top_d)[order(x, decreasing = TRUE)] }))
  d_weights <- t(apply(top_d, 1, function(x){ x[order(x, decreasing = TRUE)] }))
  d_ranks <- RankAggreg(x = d_lists,
                        seed = 0,
                        weights = d_weights,
                        k = min(ncol(d_lists), 20),
                        importance = imp,
                        standardizeWeights = FALSE,
                        verbose=FALSE)
  
  #make a nice data frame of the top from this method of aggregating.
  top_all <- d %>% mutate(comb = paste(Name, Info, cluster_size, sep = '_&_'))
  res <- t(sapply(d_ranks$top.list, function(x) top_all[top_all$comb == x,]))
  rownames(res) <- NULL
  res <- data.frame(Name = as.character(unlist(res[,1])),res[,-1], stringsAsFactors = FALSE) %>% select(-comb)
  return(res)
}

## Creates an ordered list of cluster size,
##  method name, and parameter choices for each individual measure.
getOverallLists <- function (d, n = NULL) {
  d <- d %>% unite('mpc', Name:cluster_size, remove = TRUE, sep='_&_' ) %>%
    column_to_rownames(var = 'mpc') %>%
    t()
  d <- apply(d, 1, function(x) {colnames(d)[order(x, decreasing = TRUE)]}) %>% 
    as.data.frame()
  rownames(d) <- paste('Rank', 1:nrow(d), sep = ' ')
  
  if(is.null(n)) n = ncol(d)
  d <- d[1:n,]
  
  methods <- as.data.frame(apply(d, 2, function (x) {
    unlist(lapply(strsplit(x, split = '_&_'), function(y){ y[1] }))
  }))
  params <- as.data.frame(apply(d, 2, function (x) {
    unlist(lapply(strsplit(x, split = '_&_'), function(y){ y[2] }))
  }))
  clusters <- as.data.frame(apply(d, 2, function (x) {
    unlist(lapply(strsplit(x, split = '_&_'), function(y){ y[3] }))
  }))
  
  return(list(cluster_sizes = clusters,
              methods_names = methods,
              parameter_choices = params))
}


tabler <- function (y) {
    y[,2] <- apply(y, 1, function(x) paste0(x[2]))
    y <- y %>% mutate(across(-c(Info,Name,cluster_size), ~ round(as.numeric(.), 2)))
  return(y)
}

run_lists <- function(name){
  load(paste0('measures_',name,'_formatted.RData'))
  getOverallLists(d=measures, n = 5)$cluster_sizes
}

run_lists("BAISE")
run_lists("POLLEN")
run_lists("DENG")
run_lists("GOOLAM")
run_lists("YAN")
run_lists("BEUTNER")

run_all <- function(name){
  load(paste0('measures_',name,'_formatted.RData'))
  results = getRA(d = measures, m = colnames(measures[,-c(1,2,3)]))
  res <- tabler(results)
  return(results)
  #latex_table <- xtable(res)
  #return(latex_table)
}

run_all("BAISE")
run_all("POLLEN")
run_all("DENG")
run_all("GOOLAM")
run_all("YAN")
run_all("BEUTNER")


## OPTIONAL
## Create the plot
create_plot <- function(name) {
  load(paste0('measures_', name, '_formatted.RData'))
  
  data_long <- measures %>%
    pivot_longer(cols = -c(1:3), names_to = "Measure", values_to = "Value")
  
  num_measures <- length(unique(data_long$Measure))  # Count number of facets
  
  # Dynamically set plot width/height (adjust as needed)
  plot_width <- max(7, min(14, num_measures * 2))  # 2 inches per facet, but within a reasonable range
  plot_height <- max(5, min(10, ceiling(num_measures / 3) * 3))  # 3 inches per row
  
  plot <- ggplot(data_long, aes(x = cluster_size, y = Value, shape = Name)) +
    geom_point(size = 3) +
    facet_wrap(~Measure) +
    theme_minimal(base_size = 16) +  # Increase base font size for all text elements
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      strip.text = element_text(size = 18),  # Facet labels
      axis.title = element_text(size = 18),  # Axis titles
      axis.text = element_text(size = 16),  # Axis text
      legend.title = element_text(size = 20),  # Legend title
      legend.text = element_text(size = 18)  # Legend labels
    ) +
    labs(
      title = paste0("Scaled Metric values by Cluster Size"),
      x = "Cluster Size",
      y = "Standardized Value",
      shape = "Name"
    ) +
    scale_shape_manual(values = 1:length(unique(data_long$Name)))
  
  # Save the plot
  output_path <- paste0("plot_", name, ".jpg")
  ggsave(output_path, plot, width = plot_width, height = plot_height, dpi = 300)
  
  message("Plot saved to: ", output_path)
}


create_plot("BAISE")
create_plot("POLLEN")
create_plot("DENG")
create_plot("GOOLAM")
create_plot("YAN")
create_plot("BEUTNER")

