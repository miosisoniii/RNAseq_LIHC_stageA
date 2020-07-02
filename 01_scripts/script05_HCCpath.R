#-------------------------------------------------------------------------------------#
# Project: HCC pathways - Stage A
# Purpose: Identify HCC-associated genes from literature in BG output for Stage A
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
# + Load HCC genes 
# + Load BG output for transcripts
#-------------------------------------------------------------------------------------#
hcc_genes <- read.csv("./04_HCCpath/HCCpath_genes.csv")
inputpaths <- dir("./00_data/stage1", full.names = TRUE)
race_dfs <- lapply(inputpaths, function(x) read.delim(x))

races <- c("asian", "black", "white")

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
# Filtration for HCC pathway genes through TRANSCRIPT
# + filter for hcc genes
# + retain all transcripts for 
# + write to csv
#-------------------------------------------------------------------------------------#
clean_dfs <- lapply(race_dfs, clean_deg)

allHCC_allraces <- lapply(clean_dfs, function(x) filter(x, gene_name %in% hcc_genes$name))

names(allHCC_allraces) <- races[1:3]
lapply(1:length(allHCC_allraces), function(i) write.csv(allHCC_allraces[[i]],
                                                          file = paste0("./02_output/04_HCCpathOutput/", 
                                                                        names(allHCC_allraces[i]), "_HCCretained_stage1.csv"), row.names = FALSE))

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

joinfunction(allHCC_allraces) -> allracejoined

write.csv(allracejoined, paste0("./02_output/04_HCCpathOutput/allrace_mirnabiomarkertargets_stage1.csv") )

# can include option/argument in function here to remove redundant transcripts and keep only those with Qval < 0.05 or greatest absolute F