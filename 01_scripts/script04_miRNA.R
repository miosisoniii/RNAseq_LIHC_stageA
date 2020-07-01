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
# Section: Adjust Dataset Names
#-------------------------------------------------------------------------------------#
dataname <- "stage1"
#dataname <- "stageA"

#-------------------------------------------------------------------------------------#
# Input files
# + Load PRE-filtered FC/qval cutoff bg output 
# + Load miRNA biomarkers
#-------------------------------------------------------------------------------------#
inputpaths <- dir("./00_data/stage1", full.names = TRUE)
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
                                                       file = paste0("./02_output/03_biomarkerOutput/stage1/", names(allmiRNA_allraces[i]), "_miRNA_retained_", dataname, ".csv"), row.names = FALSE))

topfc_mirna <- lapply(allmiRNA_allraces, function(x) retain_topfc(x))
inner_join(topfc_mirna [[1]],topfc_mirna [[2]], by = "gene_name", suffix = c(".asian", ".black")) -> asianblackmerge
inner_join(asianblackmerge, topfc_mirna [[3]], by = "gene_name") -> allmerge

allmerge %>% rename(pval.white = pval) %>%
  rename(abs_fc.white = abs_fc) %>%
  rename(fc.white = fc) %>%
  rename(logfc.white = logfc) %>%
  rename(qval.white =  qval) -> allmerge.final

write.csv(paste0(allmerge.final, "./02_output/03_biomarkerOutput/stage1/allrace_miRNA_", dataname, ".csv"))

#-------------------------------------------------------------------------------------#
# miRNA Biomarker targets
# Found using Xu et al. 2018, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6171016/
# + filter pre-filtered for miRNA targets
# + join lists by race
#-------------------------------------------------------------------------------------#
mirna$targets -> targets.mirna

allrace.targets <- lapply(clean_dfs, function(x) filter(x, gene_name %in% targets.mirna))

# function that takes rows with unique gene names and the lowest qval
joinfunction <- function(listobject){
  lapply(listobject, function(x) x %>% group_by(gene_name) %>% arrange(-abs_fc) %>% slice(1)) -> listobject
  inner_join(listobject[[1]], listobject[[2]], by="gene_name", suffix=c(".asian",".black")) -> ABmerge
  inner_join(ABmerge, listobject[[3]], by="gene_name") -> ABWmerge
  ABWmerge %>% rename(pval.white = pval) %>%
    rename(abs_fc.white = abs_fc) %>%
    rename(fc.white = fc) %>%
    rename(logfc.white = logfc) %>%
    rename(qval.white =  qval) -> ABWmerge
  #ABWmerge %>% group_by(gene_name)
  print(paste0("# genes for Asian: ", nrow(listobject[[1]])))
  print(paste0("# genes for Black: ", nrow(listobject[[2]])))
  print(paste0("# genes for White: ", nrow(listobject[[3]])))
  print(paste0("Total number of overlapped genes: ", nrow(ABWmerge)))
  return(ABWmerge)
}

joinfunction(allrace.targets) -> allracejoined

write.csv(allracejoined, paste0("./02_output/03_biomarkerOutput/stage1/stage1allrace_mirnabiomarkertargets.csv") )

# can include option/argument in function here to remove redundant transcripts and keep only those with Qval < 0.05 or greatest absolute F