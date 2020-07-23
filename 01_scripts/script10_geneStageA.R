#-------------------------------------------------------------------------------------#
# Project: Stage A 
# Purpose: Determine significance of Stage A gene in Treatment and Biomarkers code
# Author: Artemio Sison III
# R Version: 4.0.1 "See Things Now"
#-------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------#
# Section: Load Dependencies
#-------------------------------------------------------------------------------------#
# Install dependencies
library(dplyr)
library(stringr)

#-------------------------------------------------------------------------------------#
# Section: Read in Gene data for Stage A
#-------------------------------------------------------------------------------------#
# Read in ballgown files
inputpaths <- dir("./09_geneStageA", full.names = TRUE)
race_dfs <- lapply(inputpaths, function(x) read.delim(x))

#-------------------------------------------------------------------------------------#
# Section: Clean Data
# + Functions
# + Execute functions
#-------------------------------------------------------------------------------------#

# Select variables gene_name, fc, qval, pval
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

selectedvars <-  lapply(race_dfs, selectvar)
filtered_df <- lapply(selectedvars, filter_deg)

#-------------------------------------------------------------------------------------#
# Section: Biomarkers
# + 
#-------------------------------------------------------------------------------------#

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

biomarkerfilt <- lapply(filtered_df, function(x) filter(x, gene_name %in% biomarkers))

# create list of names for df's in biomarker_allraces object
intermediatenames <- c("asian", "black", "white")
# set names for biomarker_allraces object
names(biomarkerfilt) <- intermediatenames[1:3]
# write biomarker_allraces object to csv, to vet for transcript copies with favorable FC
lapply(1:length(biomarkerfilt), function(i) write.csv(biomarkerfilt[[i]],
                                                           file = paste0("./02_output/08_geneStageA_Output/", names(biomarkerfilt[i]), "_biomarkerRetained.csv"), row.names = FALSE))


#-------------------------------------------------------------------------------------#
# Section: Treatment
# + List treatment genes
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

lapply(filtered_df, function(x) filter(x, gene_name %in% treatment.gene)) -> treat.all


names(treat.all) <- intermediatenames[1:3]
lapply(1:length(treat.all), function(i) write.csv(treat.all[[i]],
                                                  file = paste0("./02_output/08_geneStageA_Output/", names(treat.all[i]), "_HCCtreat_retained_stageA.csv"), row.names = FALSE))
