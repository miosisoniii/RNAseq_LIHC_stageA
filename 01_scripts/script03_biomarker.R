#-------------------------------------------------------------------------------------#
# Project: Ballgown overlap
# Purpose: Identify overlapping genes stratified across race
# Author: Artemio Sison III
# R Version: 4.0.1 "Holding the Windsock"
#-------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------#
# Install dependencies
#-------------------------------------------------------------------------------------#
library(dplyr)
library(stringr)

#-------------------------------------------------------------------------------------#
# Input files
# + Load PRE-filtered FC/qval cutoff bg output 
#-------------------------------------------------------------------------------------#
inputpaths <- dir("./00_data", full.names = TRUE)
race_dfs <- lapply(inputpaths, function(x) read.delim(x))

#-------------------------------------------------------------------------------------#
# Section: Functions
# Description: Describe functions
# ################################################################################### #
# + select variables from original data set function
# + filter gene/transcript list function for MSTRG
# + retain top absolute fold change
#-------------------------------------------------------------------------------------# 
selectvar <- function(df){
  df %>% select(gene_name, fc, qval, pval)
}

filter_deg <- function(df){
  df %>% mutate(logfc = log2(fc)) %>% 
    mutate(abs_fc = abs(logfc)) %>% 
    #filter(qval < 0.05 & abs_fc > cutoff) %>% 
    arrange(desc(qval), abs_fc) %>%
    filter(str_detect(gene_name, "MSTRG*", negate = TRUE))
}

retain_topfc <- function(df){
  df %>% group_by(gene_name) %>%
    filter(abs_fc == max(abs_fc))
}

#-------------------------------------------------------------------------------------#
# Section: Biomarker filtration
# Description: filter transcript list for HCC-related biomarkers
# ################################################################################### #
# + select variables from original data set function
# + describe cutoff
# + select variables
# + filter cutoff
# + filter race gene lists by biomarkers using lapply()
# + name df's in biomarker_allrace object
# + write biomarker list with all transcript copies to csv
# + filter for transcript with greatest absolute fold change
#-------------------------------------------------------------------------------------#

#fc.cutoff <- 2.0
#fc.cutoff <- 1.0
selectedvars <- lapply(race_dfs, selectvar)
filtered_deg <- lapply(selectedvars, filter_deg)

biomarkers <- c("AFP",      # Alpha-fetoprotein https://www.genecards.org/cgi-bin/carddisp.pl?gene=AFP
                "SPP1",     # Osteopontin/Secreted phosphoprotein1 https://www.genecards.org/cgi-bin/carddisp.pl?gene=SPP1&keywords=spp1
                "MDK",      # Midkine https://www.genecards.org/cgi-bin/carddisp.pl?gene=MDK
                # GALAD Score not found, not gene
                # AFP-L3 -Alpha fetoprotein L3 not found, under AFP
                "F2",       # Prothrombin https://www.genecards.org/cgi-bin/carddisp.pl?gene=F2
                "DCP",      # Des-carboxyprothrombin https://www.gastrojournal.org/article/S0016-5085(16)31781-4/fulltext
                "DCPS",     # Decapping Enzyme, scavenger https://www.genecards.org/cgi-bin/carddisp.pl?gene=DCPS&keywords=dcp
                "DCP2",     # Des-carboxyprothrombin/decapping mRNA2 https://www.genecards.org/cgi-bin/carddisp.pl?gene=DCP2
                "DKK1",     # Dickkopf-1 https://www.genecards.org/cgi-bin/carddisp.pl?gene=DKK1
                "GPC3",     # Glypican 3 https://www.genecards.org/cgi-bin/carddisp.pl?gene=GPC3&keywords=gpc3
                "FUC1",     # alpha-1-fucosidase precursor? https://www.uniprot.org/uniprot/Q8GW72
                "FUCA1",    # alpha-L-fucosidase 1 https://www.genecards.org/cgi-bin/carddisp.pl?gene=FUCA1
                "GOLM1",    # golgi protein 1 https://www.uniprot.org/uniprot/Q8NBJ4
                "GP73",     # golgi protein 73, golph2 https://onlinelibrary.wiley.com/doi/full/10.1111/cpr.12538
                "SCCA",     # squamous cell carcinoma antigen https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3079753/
                "SCCA1",    # squamous cell carcinoma antigen variant 1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3079753/
                "SCCA2")    # squamous cell carcinoma antigen variant 2 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3079753/


biomarker_allraces <- lapply(filtered_deg, function(x) filter(x, gene_name %in% biomarkers))

# create list of names for df's in biomarker_allraces object
intermediatenames <- c("asian", "black", "white")
# set names for biomarker_allraces object
names(biomarker_allraces) <- intermediatenames[1:3]
# write biomarker_allraces object to cvs, to vet for transcript copies with favorable FC
lapply(1:length(biomarker_allraces), 
       function(i) write.csv(biomarker_allraces[[i]],
                             file = paste0("./02_output/03_biomarkerOutput/", 
                                           names(biomarker_allraces[i]), "_biomarkerRetained.csv"), row.names = FALSE))


topfc_allraces <- lapply(biomarker_allraces, function(x) retain_topfc(x))


inner_join(biomarker_topfc[[1]],biomarker_topfc[[2]], by = "gene_name", suffix = c(".asian", ".black")) -> asianblackmerge
left_join(asianblackmerge, biomarker_topfc[[3]], by = "gene_name") -> allmerge

allmerge %>% rename(pval.white = pval) %>%
  rename(abs_fc.white = abs_fc) %>%
  rename(fc.white = fc) %>%
  rename(logfc.white = logfc) %>%
  rename(qval.white =  qval) -> allmerge.final



write.csv(results_genes, paste0("./", race, "/02_ballgown/resulttables/02_genenames.csv"))
