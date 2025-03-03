#===================================================================================================
# 01_format_data.R is meant to format each single cell dataset into the same format:
# - rows with genes
# - columns with samples
# - entries with raw counts
# The formatted matrices for each dataset are then saved to another directory
#===================================================================================================

# Packages
library(org.Mm.eg.db)
library(dplyr)
library(tibble)
library(purrr)
library(data.table)
library(R.utils)

# manually set session to data file path

# Here is the index
# as.data.frame(list.files())
# > wd.list
# 1             Baise_sc_data
# 2  Beutner_mice_age_sc_data
# 3              Deng_sc_data
# 4            Goolam_sc_data
# 8            Pollen_sc_data
# 13              Yan_sc_data

## BAISE_sc_data
# Read the tsv into tables
BAISE1 = read.table(file = "Baise_sc_data/FPKM-2cell-upper-quartile-normalized.txt", header=T, row.names = 1, sep = "\t")
BAISE2 = read.table(file = "Baise_sc_data/FPKM-4cell-upper-quartile-normalized.txt", header=T, row.names = 1, sep = "\t") 
BAISE3 = read.table(file = "Baise_sc_data/FPKM-zygotes-upper-quartile-normalized.txt", header=T, row.names = 1, sep = "\t")
BAISE.meta = c(rep('2cell', length(BAISE1)),rep('4cell', length(BAISE2)),rep('zygote', length(BAISE3)))

# remove pointless names
BAISE1 = BAISE1[-1,-1] %>% rownames_to_column(var = "Genes")
BAISE2 = BAISE2[-1,-1] %>% rownames_to_column(var = "Genes")
BAISE3 = BAISE3 %>% rownames_to_column(var = "Genes")


# merge cells together
BAISE_data = Reduce(function(x, y) merge(x, y, by="Genes", all = TRUE), list(BAISE1, BAISE2, BAISE3)) %>%
  column_to_rownames(var = "Genes") 
colnames(BAISE_data) <- sapply(1:ncol(BAISE_data), function(x) paste0(BAISE.meta[x], '_cell', x))
save(BAISE_data, file = "Baise_sc_data/BAISE_data.RData")


##  Deng_sc_data
# RUN ONCE
# for(f in list.files("Deng_sc_data/extracted_files")){
#   for(g in list.files(paste0("Deng_sc_data/extracted_files/", f))){
#     gunzip(paste0("Deng_sc_data/extracted_files/", f,"/",g), remove = TRUE) }}

# Go through the unzipped files and get genes and reads
cells <- list()
for(f in list.files("Deng_sc_data/extracted_files")){
  for(g in list.files(paste0("Deng_sc_data/extracted_files/", f))){
    cell_data <- read.table(file = paste0("Deng_sc_data/extracted_files/", f,"/",g),
                            sep="\t")
    cells[[g]] <- cell_data[,c(1,4)] }}

# Bind all the cells together and set row/column names
DENG_data <- do.call(cbind, lapply(cells, function(x) x[,2]))
meta.info <- c(rep('16cell', 50),  rep('4cell', 14),  rep('8cell', 37),  rep('early2cell', 8),  rep('earlyBlast', 43),
               rep('late2cell', 10),  rep('lateBlast',30),  rep('mid2cell',12),  rep('midBlast',60),  rep('zygote',4))
colnames(DENG_data) <- sapply(1:ncol(DENG_data), function(x) paste0(meta.info[x], '_cell', x))
rownames(DENG_data) <- cells[[1]][,1]
save(DENG_data, file = "Deng_sc_data/DENG_data.RData")


## Goolam_sc_data
list.files("Goolam_sc_data")
GOOLAM_data <- read.table(file = "Goolam_sc_data/Goolam_et_al_2015_count_table.tsv", header=T, sep = "\t")
GOOLAM.meta <- colnames(GOOLAM_data)
c2 <- which(grepl('X2cell', GOOLAM.meta)); GOOLAM.meta[c2] <- 'X2cell'
c4 <- which(grepl('4cell', GOOLAM.meta)); GOOLAM.meta[c4] <- '4cell'
c8 <- which(grepl('8cell', GOOLAM.meta)); GOOLAM.meta[c8] <- '8cell'
c16 <- which(grepl('16cell', GOOLAM.meta)); GOOLAM.meta[c16] <- '16cell'
c32 <- which(grepl('32cell', GOOLAM.meta)); GOOLAM.meta[c32] <- '32cell'
colnames(GOOLAM_data) <- sapply(1:ncol(GOOLAM_data), function(x) paste0(GOOLAM.meta[x], '_cell', x))
save(GOOLAM_data, file = "Goolam_sc_data/GOOLAM_data.RData")


## Pollen_sc_data
list.files("Pollen_sc_data")
POLLEN_data <- read.table(file = "Pollen_sc_data/HiSeq301-RSEM-linear values.txt", header=T, sep = "\t")
POLLEN_data <- POLLEN_data %>% column_to_rownames(var = "Gene_Symbol")
POLLEN.meta = sapply(colnames(POLLEN_data), function(x) sub("\\_[^_]*$", "", x) )
colnames(POLLEN_data) <- sapply(1:ncol(POLLEN_data), function(x) paste0(POLLEN.meta[x], '_cell', x))
save(POLLEN_data, file = "Pollen_sc_data/POLLEN_data.RData")


## Yan_sc_data
list.files("Yan_sc_data")
YAN_data <- readxl::read_xlsx("Yan_sc_data/yanstuff.xlsx")
YAN_data <- YAN_data[,-2] %>% column_to_rownames(var="Gene_ID")

YAN.meta <- colnames(YAN_data)
YAN.meta[which(grepl('Oocyte', YAN.meta))] <- 'Oocyte'
YAN.meta[which(grepl('Zygote', YAN.meta))] <- 'Zygote'
YAN.meta[which(grepl('2-cell', YAN.meta))] <- '2cell'
YAN.meta[which(grepl('4-cell', YAN.meta))] <- '4cell'
YAN.meta[which(grepl('8-cell', YAN.meta))] <- '8cell'
YAN.meta[c(which(grepl('blast', YAN.meta)), which(YAN.meta == 'Morulae #1 -Cell#3(RPKM)'), which(YAN.meta == 'Morulae #1 -Cell#8(RPKM)') ) ] <- 'Blastocyst'
YAN.meta[which(grepl('Morulae', YAN.meta))] <- 'Morula'
YAN.meta[which(grepl('passage', YAN.meta))] <- 'hESCpassage'
colnames(YAN_data) <- sapply(1:ncol(YAN_data), function(x) paste0(YAN.meta[x], '_cell', x))
save(YAN_data, file = "Yan_sc_data/YAN_data.RData")


## Beutner_mice_age_sc_data
list.files("Beutner_mice_age_sc_data")
gunzip("Beutner_mice_age_sc_data/GSM4318801_3M_matrix.txt.gz", remove = TRUE)
BEUTNER_data <- read.table("Beutner_mice_age_sc_data/GSM4318801_3M_matrix.txt", header=TRUE)

# get gene info
gene_info <- select(org.Mm.eg.db, keys = rownames(BEUTNER_data), keytype = 'ALIAS', columns = c('ALIAS', "GENENAME"))
gene_info <- cbind(gene_info, index = 1:nrow(gene_info))
gene_info <- gene_info %>% group_by(ALIAS) %>% summarise(index = min(index), GENENAME  = paste(GENENAME, collapse = ", ")) %>% arrange(index)

head(gene_info)

# find mitochondrial cells so we can filter low quality
mitoch_ind <- grepl("mitoch", gene_info$GENENAME)
sum(mitoch_ind)

# calculate %mt of all cells. 
cell_mt_percent <- apply(BEUTNER_data, 2, function(x) {  sum(x[mitoch_ind]) / sum(x)  })
cell_gene_count <- apply(BEUTNER_data, 2, sum)

# doublets or cells with poor quality (genes > 6000, genes < 200, or >5% genes mapping to mitochondrial genome) are excluded
BEUTNER_data <- BEUTNER_data[, cell_mt_percent < 0.05
                                 & 
                                   cell_gene_count > 200
                                 & 
                                   cell_gene_count < 6000]

# this leaves us about 2k cells.
save(BEUTNER_data, file = "Beutner_mice_age_sc_data/BEUTNER_data.RData")

