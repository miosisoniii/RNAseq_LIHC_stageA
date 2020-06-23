#-------------------------------------------------------------------------------------#
# Project: miRNA identification - Stage A
# Purpose: Identify miRNA from literature in BG output for Stage A
# Author: Artemio Sison III
# R Version: 4.0.1 "Holding the Windsock"
#-------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------#
# Install dependencies
#-------------------------------------------------------------------------------------#
library(dplyr)
library(stringr)
library(tidyr)

#-------------------------------------------------------------------------------------#
# Input files
# + Load PRE-filtered FC/qval cutoff bg output 
# + Load miRNA biomarkers
#-------------------------------------------------------------------------------------#
inputpaths <- dir("./00_data", full.names = TRUE)
race_dfs <- lapply(inputpaths, function(x) read.delim(x))

races <- c("asian", "black", "white")

mirna <- read.csv("./03_miRNA/xu_miRNAbiomarkers.csv")
mirna <- mirna[1:3]

#-------------------------------------------------------------------------------------#
# Section: Functions
# Description: Describe functions
# ################################################################################### #
# + clean BG output
# + retain top absolute fold change
#-------------------------------------------------------------------------------------# 
clean_deg <- function(df){
  df %>% select(gene_name, fc, qval, pval) %>%
    mutate(logfc = log2(fc)) %>% mutate(abs_fc = abs(logfc)) %>% 
    #filter(qval < 0.05 & abs_fc > cutoff) %>% # ignore cutoff for miRNA identification
    arrange(desc(qval), abs_fc) %>%
    filter(str_detect(gene_name, "MSTRG*", negate = TRUE))
}

retain_topfc <- function(df){
  df %>% group_by(gene_name) %>%
    filter(abs_fc == max(abs_fc))
}

#-------------------------------------------------------------------------------------#
# Data Cleanup -miRNA
# Using https://tidyr.tidyverse.org/reference/unnest.html
# + un-nest multiple miRNA
# + function to create logfc, abs_fc, qval, pval
# + function to select appropriate variables
#-------------------------------------------------------------------------------------#

mirna$biomarker <- as.character(mirna$biomarker)
mirna$targets <- as.character(mirna$targets)
mirna %>% unnest(targets = strsplit(targets, ","), keep_empty = TRUE) %>% unnest(biomarker = strsplit(biomarker, ","), keep_empty = TRUE)-> mirna

str_replace(mirna$biomarker, "miR-", "MIR") -> mirna$biomarker
str_replace(mirna$targets, "No mentioned", "NA") -> mirna$targets

write.csv(mirna, "./03_miRNA/mirna_clean.csv")

#-------------------------------------------------------------------------------------#
# Filter for miRNA
# Found using Xu et al. 2018, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6171016/
# + filter pre-filtered for miRNA targets NOTE: NONE FOUND IN ANY RACE
# + join lists by race
# + identify miRNA in datasets with str_detect()
#-------------------------------------------------------------------------------------#

# list of transcripts do NOT contain any miRNA
#mirna$biomarker -> mirna_biomarker

clean_dfs <- lapply(race_dfs, clean_deg)

allmiRNA_allraces <- lapply(clean_dfs, function(x) filter(x, str_detect(gene_name, "MIR")))

names(allmiRNA_allraces) <- races[1:3]
lapply(1:length(allmiRNA_allraces), function(i) write.csv(allmiRNA_allraces[[i]],
                                                       file = paste0("./03_miRNA/", names(allmiRNA_allraces[i]), "_miRNA_retained.csv"), row.names = FALSE))

topfc_mirna <- lapply(allmiRNA_allraces, function(x) retain_topfc(x))
inner_join(topfc_mirna [[1]],topfc_mirna [[2]], by = "gene_name", suffix = c(".asian", ".black")) -> asianblackmerge
inner_join(asianblackmerge, topfc_mirna [[3]], by = "gene_name") -> allmerge

allmerge %>% rename(pval.white = pval) %>%
  rename(abs_fc.white = abs_fc) %>%
  rename(fc.white = fc) %>%
  rename(logfc.white = logfc) %>%
  rename(qval.white =  qval) -> allmerge.final

write.csv(allmerge.final, "./03_miRNA/allrace_miRNA.csv")
