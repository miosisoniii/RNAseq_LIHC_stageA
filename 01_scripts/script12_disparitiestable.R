#-------------------------------------------------------------------------------------#
# Project: Table Generation - LIHC Stage 1 Cancer Disparities
# Purpose: Create table for describing disparities in race
# Author: Artemio Sison III
# R Version: 4.0.1 "See Things Now"
#-------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------#
# Load Depencies
#-------------------------------------------------------------------------------------#
rm(list = ls())
library(dplyr)

#-------------------------------------------------------------------------------------#
# Load DEG data
#-------------------------------------------------------------------------------------#

# Read in ballgown files
inputpaths <- dir("./00_data/stage1", full.names = TRUE)
race_dfs <- lapply(inputpaths, function(x) read.delim(x))

#-------------------------------------------------------------------------------------#
# Functions
#-------------------------------------------------------------------------------------#

# Filter for DEG with cutoff and create logFC and absolute FC, select highest absolute FC
filter_deg <- function(df, cutoff){
  df %>% mutate(logfc = log2(fc)) %>% 
    mutate(abs_fc = abs(logfc)) %>% 
    filter(qval < 0.05 & abs_fc > cutoff) %>% 
    arrange(desc(qval), abs_fc) %>%
    filter(str_detect(gene_name, "MSTRG*", negate = TRUE))
}

# retain lowest Q value (or highest absolute fold change)
retain_singleDEG <- function(df){
  df %>% group_by(gene_name) %>%
    filter(qval == min(qval))
}

# join race gene tables together
joinfunction <- function(listobject){
  lapply(listobject, function(x) x %>% group_by(gene_name) %>% arrange(-abs_fc) %>% slice(1)) -> listobject
  inner_join(listobject[[1]], listobject[[2]], by="gene_name", suffix=c(".asian",".black")) -> ABmerge
  inner_join(ABmerge, listobject[[3]], by="gene_name") -> ABWmerge
  ABWmerge %>% 
    rename(pval.white = pval) %>%
    rename(abs_fc.white = abs_fc) %>%
    rename(fc.white = fc) %>%
    rename(logfc.white = logfc) %>%
    rename(qval.white =  qval) -> ABWmerge
  print(paste0("# genes for Asian: ", nrow(listobject[[1]])))
  print(paste0("# genes for Black: ", nrow(listobject[[2]])))
  print(paste0("# genes for White: ", nrow(listobject[[3]])))
  print(paste0("Total number of overlapped genes: ", nrow(ABWmerge)))
  return(ABWmerge)
}

# create opposite of %in% to select genes not inside gene list
"%ni%" <- Negate("%in%")

#-------------------------------------------------------------------------------------#
# Clean data
#-------------------------------------------------------------------------------------#

# Set cutoff for logFC
cutoff <- 1.0

# access dataframe list and apply filter
filtered_deg <- lapply(race_dfs, filter_deg, cutoff)

# retain single DEG
filtered_deg <- lapply(filtered_deg, retain_singleDEG)

#-------------------------------------------------------------------------------------#
# Identify genes unique to each race
#-------------------------------------------------------------------------------------#

# merge all dfs into one df
overlapgenes <- joinfunction(filtered_deg)

# select overlapping genes to remove from anti join comparisons for each race
overlapgenelist <- overlapgenes$gene_name

# remove overlapped genes from combined genes list for all 3 races
overlapremoved <- lapply(filtered_deg, function(x) x %>% filter(gene_name %ni% overlapgenelist))

## use %ni% to see if there are unique genes present in asian and not white
# assign list of genes
asiangenelist <- overlapremoved[[1]]$gene_name
whitegenelist <- overlapremoved[[3]]$gene_name

# use list to "anti-filter" for genes not common in respective list
unique_asian <- overlapremoved[[1]] %>% filter(gene_name %ni% whitegenelist)
unique_white <- overlapremoved[[3]] %>% filter(gene_name %ni% asiangenelist)


#-------------------------------------------------------------------------------------#
# Create unique Upregulated DEG table
#-------------------------------------------------------------------------------------#

# identify upregulated genes for asian
upasian_unique <- unique_asian %>% select(gene_name, logfc, qval, pval) %>% 
  mutate(race = "asian") %>%
  mutate(exp = "up") %>%
  filter(logfc > 0) %>% arrange(qval)

# identify upregulated genes for black 
upblack_unique <- overlapremoved[[2]] %>% select(gene_name, logfc, qval, pval) %>%
  mutate(race = "black") %>%
  mutate(exp = "up") %>%
  filter(logfc > 0) %>% arrange(qval)

# identify upregulated genes for white
upwhite_unique <- unique_white %>% select(gene_name, logfc, qval, pval) %>% 
  mutate(race = "white") %>%
  mutate(exp = "up") %>%
  filter(logfc > 0) %>% arrange(qval)

#-------------------------------------------------------------------------------------#
# Create unique Downregulated DEG table
#-------------------------------------------------------------------------------------#

# identify downregulated genes for asian
downasian_unique <- unique_asian %>%  select(gene_name, logfc, qval, pval) %>% 
  mutate(race = "asian") %>%
  mutate(exp = "down") %>%
  filter(logfc < 0) %>% arrange(qval)

# identify downregulated genes for black 
downblack_unique <- overlapremoved[[2]] %>%  select(gene_name, logfc, qval, pval) %>% 
  mutate(race = "black") %>%
  mutate(exp = "down") %>%
  filter(logfc < 0) %>% arrange(qval)

# identify downregulated genes for white
downwhite_unique <- unique_white %>%  select(gene_name, logfc, qval, pval) %>%
  mutate(race = "white") %>%
  mutate(exp = "down") %>%
  filter(logfc < 0) %>% arrange(qval)

#-------------------------------------------------------------------------------------#
# Select top 10 genes based on qvalue, combine datasets
#-------------------------------------------------------------------------------------#

# upregulated gene lists
# take top 10 genes for each race
top_upasian <- head(upasian_unique, n = 10) 
top_upblack <- head(upblack_unique, n = 10)
top_upwhite <- head(upwhite_unique, n = 10)

# rowbind into same dataset
upreg <- bind_rows(top_upasian, top_upblack, top_upwhite)

# downregulated gene lists
# take top 10 genes for each race
top_downasian <- head(downasian_unique, n = 10) 
top_downblack <- head(downblack_unique, n = 10)
top_downwhite <- head(downwhite_unique, n = 10)

# rowbind into same dataset
downreg <- bind_rows(top_downasian, top_downblack, top_downwhite)

#-------------------------------------------------------------------------------------#
# write to csv
#-------------------------------------------------------------------------------------#
write.csv(upreg, "./02_output/10_diffgenes/upreg_fc1.csv")
write.csv(downreg, "./02_output/10_diffgenes/downregfc2.csv")
