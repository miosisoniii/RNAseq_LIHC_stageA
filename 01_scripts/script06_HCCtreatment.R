#-------------------------------------------------------------------------------------#
# Project: HCC treatement identification - Stage 1
# Purpose: Identify HCC treatment genes from literature in BG output for Stage 1
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
# + Load BG output for transcripts
#-------------------------------------------------------------------------------------#
inputpaths <- dir("./00_data/stage1", full.names = TRUE)
race_dfs <- lapply(inputpaths, function(x) read.delim(x))
races <- c("asian", "black", "white")


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
# Section: Match differentially expressed genes to treatment
#-------------------------------------------------------------------------------------#
treatment.gene <- c("ARAF",      # RAF KINASES (ARAF, BRAF, CRAF)
                    "RAF1",      # RAF KINASE?
                    "FLT1",      # VEGFR1 aka FLT1?
                    "KDR",       # VEGFR2 aka KDR?
                    # VEGFR3 aka FLT4
                    "VEGFB",     # VEGF - FOUND VEGFB, VEGFC
                    "VEGFC",     # VEGF - FOUND VEGFB, VEGFC
                    "PDGFRB",    # PDGFRB
                    "PDGFRA",    # PDGFRA
                    "FGFR1OP2",  # FGFR1 - FOUND FGFR1OP2 https://www.genecards.org/cgi-bin/carddisp.pl?gene=FGFR1OP2
                    "FGFR2",     # FGFR2
                    "FGFR3",     # FGFR3
                    "FGFR4",     # FGFR4
                    "FGFRL1",    # FGFR1 - FOUND FGFRL1?
                    "RET",       # RET
                    "TIE1",      # TIE2 - FOUND TIE1?
                    "TEKT4P2",   # TIE2 aka TEK - FOUND TEKT4P2
                    # NOT FOUND  # KIT
                    "PDCD10",    # PD-1 aka PDCD1 - FOUND PDCD10
                    "PDCD11",    # PD-1 aka PDCD1 - FOUND PDCD11
                    "PDCD4",     # PD-1 aka PDCD1 - FOUND PDCD4
                    "PDCD5P1",   # PD-1 aka PDCD1 - FOUND PDCD5P1 https://www.genecards.org/cgi-bin/carddisp.pl?gene=PDCD5P1
                    "PDCD7",     # PD-1 aka PDCD1 - FOUND PDCD7
                    "PDCD5",     # PD-1 aka PDCD1 - FOUND PDCD5
                    "PDCD2L",    # PD-1 aka PDCD1 - FOUND PDCD2L
                    "PDCD6IP",   # PD-1 aka PDCD1 - FOUND PDCD6IP
                    "PDCD6",     # PD-1 aka PDCD1 - FOUND PDCD6
                    "PDCD2L",    # PD-1 aka PDCD1 - FOUND PDCD2L
                    "MET"        # C-Met aka MET
)    

race.all <- lapply(race_dfs, clean_deg)

lapply(race.all, function(x) filter(x, gene_name %in% treatment.gene)) -> treat.all


names(treat.all) <- races[1:3]
lapply(1:length(treat.all), function(i) write.csv(treat.all[[i]],
                                                          file = paste0("./02_output/05_HCCtreatOutput/", names(treat.all[i]), "_HCCtreat_retained_stage1.csv"), row.names = FALSE))


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

joinfunction(treat.all) -> allracejoined

write.csv(allracejoined, paste0("./02_output/05_HCCtreatOutput/allrace_HCCtreat_stage1.csv") )

    