#-------------------------------------------------------------------------------------#
# Project: Ballgown overlap
# Purpose: Identify overlapping genes stratified across race
# Author: Artemio Sison III
# R Version: 4.0.1 "Holding the Windsock"
#-------------------------------------------------------------------------------------#

# Install dependencies
library(dplyr)
library(stringr)

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
    filter(str_detect(gene_name, "MSTRG*", negate = TRUE))
}
# MUST CHECK TO SEE IF LOWER FC IS REMOVED!

#-------------------------------------------------------------------------------------#
# Section: Overlap DEG
# Description: Overlap DEG stratified by Race using functions
# ################################################################################### #
# + select appropriate variables
# + filter for DEG with FC > 2.0, create absFC, logFC, select genes with highest FC
#-------------------------------------------------------------------------------------#
fc.cutoff <- 1.0

selectedvars <- lapply(race_dfs, selectvar) 

filtered_deg <- lapply(selectedvars, filter_deg, fc.cutoff)

filtered_deg[[1]] -> asiandeg
filtered_deg[[2]] -> blackdeg
filtered_deg[[3]] -> whitedeg

inner_join(asiandeg, whitedeg, by = "gene_name", suffix = c(".asian", ".white")) -> joined
inner_join(joined, blackdeg, by = "gene_name") -> joined.all

joined.all %>%
  rename(logfc.black = logfc) %>%
  rename(abs_fc.black = abs_fc) %>%
  rename(fc.black = fc) %>%
  rename(pval.black = pval) %>%
  rename(qval.black = qval) -> renamed.all

renamed.all %>% distinct(gene_name, .keep_all = TRUE) -> renamed.all

write.csv(asiandeg, "./02_output/asian_degStageA_FC1.csv")
write.csv(blackdeg, "./02_output/black_degStageA_FC1.csv")
write.csv(whitedeg, "./02_output/white_degStageA_FC1.csv")

write.csv(renamed.all, "./02_output/allrace_overlapDEG_stageA_FC1.csv")

#FC 2.0
# write.csv(asiandeg, "./02_output/asian_degStageA_FC2.csv")
# write.csv(blackdeg, "./02_output/black_degStageA_FC2.csv")
# write.csv(whitedeg, "./02_output/white_degStageA_FC2.csv")
# 
# write.csv(renamed.all, "./02_output/allrace_overlapDEG_stageA_allrace.csv")

