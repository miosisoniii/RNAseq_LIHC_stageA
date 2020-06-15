#-------------------------------------------------------------------------------------#
# Project: Ballgown overlap
# Purpose: Identify overlapping genes stratified across race
# Author: Artemio Sison III
# R Version: 4.0.1 "Holding the Windsock"
#-------------------------------------------------------------------------------------#

# Install dependencies
library(dplyr)

# Read in ballgown files
inputpaths <- dir("./00_data", full.names = TRUE)
race_dfs <- lapply(inputpaths, function(x) read.delim(x))


# Functions
# Select variables gene_name, fc, qval, pval
selectvar <- function(df){
  df %>% select(gene_name, fc, qval, pval)
}

# Filter for DEG with cutoff and create logFC and absolute FC, select highest absolute FC
filter_deg <- function(df, cutoff){
  df %>% mutate(logfc = log2(fc)) %>% 
    mutate(abs_fc = abs(logfc)) %>% 
    filter(qval < 0.05 & abs_fc > cutoff) %>% 
    arrange(desc(qval), abs_fc) %>% 
    filter(abs_fc == max(abs_fc))
}


#-------------------------------------------------------------------------------------#
# Section: Overlap DEG
# Description: Overlap DEG stratified by Race using functions
# ################################################################################### #
# + select appropriate variables
# + filter for DEG with FC > 2.0, create absFC, logFC, select genes with highest FC
#-------------------------------------------------------------------------------------#
selectedvars <- lapply(race_dfs, selectvar) 

filtered_deg <- lapply(selectedvars, filter_deg, 2.0)

filtered_deg[[1]] -> asiandeg
filtered_deg[[2]] -> blackdeg
filtered_deg[[3]] -> whitedeg




