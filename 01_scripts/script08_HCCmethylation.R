#-------------------------------------------------------------------------------------#
# Project: Hepatocyte-related gene identification
# Purpose: Identify hepatocyte=specific genes from literature in BG output for Stage 1
# Author: Artemio Sison III
# R Version: 4.0.1 "Holding the Windsock"
#-------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------#
# Install dependencies
#-------------------------------------------------------------------------------------#
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

#-------------------------------------------------------------------------------------#
# Input files
# + Load BG output for transcripts
#-------------------------------------------------------------------------------------#
inputpaths <- dir("./00_data/stage1", full.names = TRUE)
race_dfs <- lapply(inputpaths, function(x) read.delim(x))
races <- c("asian", "black", "white")

#-------------------------------------------------------------------------------------#
# Section: Functions
# + clean differentially expressed genes
# + retain top FC differentially expressed genes
# + merge race df's function
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

#-------------------------------------------------------------------------------------#
# Section: Clean data
#-------------------------------------------------------------------------------------#
race.all <- lapply(race_dfs, clean_deg)

#-------------------------------------------------------------------------------------#
# Section: Filter for Methylation-related genes
# + Create list of genes associated with HCC-specific methylation
#-------------------------------------------------------------------------------------#

methylation_add.gene <- c("TBC1D1", # Original
                          "SLC39A12", # Original
                          "SLC39A11", # similar to SCL39A12
                          "SERHL", # Original 
                          "SERHL2",
                          "C6orf206", # Original
                          "MRPS18A", # https://www.genecards.org/cgi-bin/carddisp.pl?gene=RSPH9&keywords=c6orf206
                          "C11orf2", # Original
                          "IFT46",
                          "CCDC37", # Original
                          "LILRA1", # Original
                          "OR51B4", # Original
                          "HDAC1", # Original
                          "NHEJ1", # "XLF", # Original
                          "KCTD4", # Original
                          "INA", # Original
                          "ABHD9", # Original
                          "AQP6", # Original
                          "ANKRD33", # Original
                          "TSPYL5", # Original
                          "RPS6KC1", # Original
                          "FLJ13149", # Original
                          "FASTKD5",
                          "FCAR", # Original
                          "CD1B", # Original
                          "CD1A", # https://www.genecards.org/cgi-bin/carddisp.pl?gene=CD1B&keywords=CD1B
                          "KLK2", # Original
                          "HK2",
                          "CDH18", # Original
                          "CDH24",
                          "BOLL", # Original
                          "SGNE1", # Original https://www.genecards.org/cgi-bin/carddisp.pl?gene=SCG5&keywords=SGNE1
                          "SCG5",
                          "CYP11B1", # Original
                          "CPN1", # https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYP11B1&keywords=CYp11B1
                          "S100A8", # Original
                          "MAPKAP1", # Original
                          "QTRT1", # Original
                          "KLK9", # Original
                          "KLKL3",
                          "FLJ46481", # Original
                          "C4orf50", # https://www.genecards.org/cgi-bin/carddisp.pl?gene=C4orf50&keywords=FLJ46481
                          
                          ## Zhang paper methylation 
                          # Carcinoma tissue/Adjacent Tissue
                          "CDKN2A",# "P16", # HCC serum/normal serum
                          "RASSF1",#"RASSF1A", # HCC serum/normal serum
                          "TFDP3", # "APC", https://www.genecards.org/cgi-bin/carddisp.pl?gene=TFDP3
                          "GSTP1", # HCC serum/normal serum
                          "CDH1", # HCC serum/normal serum
                          "CDKN2B", #"p15", CDKN2B
                          "RUNX3", # HCC serum/normal serum
                          "SOCS1", 
                          "MGMT",
                          "PRDM2",
                          "RIZ",
                          "KMT8",
                          "SFRP1",
                          "DAPK1",
                          "DAPK2", # https://www.genecards.org/cgi-bin/carddisp.pl?gene=DAPK1&keywords=DAPk1
                          "DAPK3", # https://www.genecards.org/cgi-bin/carddisp.pl?gene=DAPK1&keywords=DAPk1
                          "CDKN2A", # "p14",
                          "CDK2NAIPNL", # additional CDK2NA
                          "RARB",
                          "ARRB2", # related to RARB https://www.genecards.org/cgi-bin/carddisp.pl?gene=ARRB2
                          "IGF2",
                          "IGF2R",
                          "TP73", # "P73", TP73
                          "MLH1",   # "hMLH1", MLH1
                          "G3BP2", # GAT3 Binding Protein G3B https://www.genecards.org/cgi-bin/carddisp.pl?gene=PRDM2&keywords=prdm2
                          "G3BP1", # GAT3 Binding Protien G3B https://www.genecards.org/cgi-bin/carddisp.pl?gene=PRDM2&keywords=prdm2
                          "WIF1", # HCC serum/normal serum
                          "SPINT2",
                          "RB1",
                          "OPCML",
                          "OBCAM", # https://www.genecards.org/cgi-bin/carddisp.pl?gene=OPCML&keywords=opcml
                          "OPCM",  # https://www.genecards.org/cgi-bin/carddisp.pl?gene=OPCML&keywords=opcml
                          "WT1", # https://www.genecards.org/cgi-bin/carddisp.pl?gene=WT1&keywords=wt1
                          "WT33", # https://www.genecards.org/cgi-bin/carddisp.pl?gene=WT1&keywords=wt1
                          "NPHS4", # https://www.genecards.org/cgi-bin/carddisp.pl?gene=WT1&keywords=wt1
                          "AWT1", # https://www.genecards.org/cgi-bin/carddisp.pl?gene=WT1&keywords=wt1
                          "WAGR", # https://www.genecards.org/cgi-bin/carddisp.pl?gene=WT1&keywords=wt1
                          "GUD", # https://www.genecards.org/cgi-bin/carddisp.pl?gene=WT1&keywords=wt1
                          
                          # Fan 2018 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6142709/
                          "PEG10",
                          "SPINK1",
                          "CDKN3",
                          "PTTG1",
                          "CENPF",
                          "ROBO1",
                          "ESM1",
                          "DTL",
                          "ANLN",
                          "MMP12",
                          "TK1",
                          "NEK2",
                          "CCL20",
                          "TREM2",
                          "MEP1A",
                          "SLC6A8",
                          "CAP2",
                          "NRXN3",
                          "RBM24",
                          
                          "MT1F",
                          "SLC22A1",
                          "SRD5A2",
                          "MT1X",
                          "BMPER",
                          "VIPR1",
                          "ESR1",
                          "HHIP",
                          "GPM6A",
                          "KCCN2",
                          "PCDH9",
                          "MYRIP",
                          "NGFR",
                          "TSLP",
                          
                          #Cheng 2018
                          "NEBL",
                          "FAM55C",
                          "GALNT3",
                          "DSE",
                          "GALNT3",
                          "DSE",
                          
                          # Tao 2011
                          # Methylation, hepatocyte-specific genes
                          "WNK2",
                          "FOXE3",
                          "MAGEA6",
                          "GUCY1A2",
                          "MILIN",
                          "TRIM58",
                          "DKFZp434I1020",
                          "ADCY5",
                          "GZMB",
                          "OVOL1",
                          "INS",
                          "ADAM8",
                          "SERPINB3",
                          "GRASP",
                          "ZMYND10",
                          "UTF1",
                          "BOLL",
                          "TM6SF1",
                          "MAGEA3",
                          "KIF17",
                          "SYK",
                          "EPHA4",
                          "ATP8A2",
                          "HIST1H4F",
                          "TLX3",
                          "NTRK3",
                          "PCDHGA12",
                          "VAV3",
                          "CELSR3",
                          "HOXB4",
                          "IRAK3",
                          "ARMCX4",
                          "HIST3H2BB",
                          "C6orf150",
                          "JAKMIP1",
                          "GULP1",
                          "NR2E1",
                          "KRTHB5",
                          "GEFT",
                          "EFCAB1",
                          "CST2",
                          "NPY",
                          "CDKL2",
                          "ZMYND10",
                          "CYSLTR1",
                          "DNAHL1",
                          "ADCY4",
                          "KRTHB3",
                          "ZAR1",
                          "GATA5",
                          "EGR2",
                          "CIAS1",
                          "FLJ42486",
                          "CNIH3",
                          "DEFA3",
                          "DRD4",
                          "KCNC4",
                          "OLFM2",
                          "NGFR",
                          "DAB2IP",
                          "ST8SIA2",
                          "SCN3B",
                          "DLEC1"
                          
)    

# filter for all HCC-specific methylation related genes
lapply(race.all, function(x) filter(x, gene_name %in% methylation_add.gene)) -> methylation.all

# retain and write to csv
names(methylation.all) <- races[1:3]
lapply(1:length(methylation.all), function(i) write.csv(methylation.all[[i]],
                                                        file = paste0("./02_output/06_SChepatocyteOutput/", names(methylation.all[i]), "_HCCmethylation_retained_stage1.csv"), row.names = FALSE))


joinfunction(methylation.all) -> HCCmeth.all