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

SCdata <- read.csv("./06_hepatocyte/scranseq_genelist.csv")
SCdata %>% rename(gene_name = gene) -> SCdata

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
# Section: filter for hepatocyte genes
#-------------------------------------------------------------------------------------#
# Tao 2011
# Methylation, hepatocyte-specific genes
hepatocytegenes <- c("WNK2",
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
                         "DLEC1")

# create data frame
hepatocytegenedf <- as.data.frame(hepatocytegenes)

# rename column to gene_name
names(hepatocytegenedf)[1] <- c("gene_name")

# row bind Tao genes (hepatocytegenedf) to SC data  df using rbind.fill from plyr
plyr::rbind.fill(SCdata, hepatocytegenedf) -> SCdata

#-------------------------------------------------------------------------------------#
# Section: Merge Single-Cell Hepatocyte-specific data with race gene data object (1/2)
# + NOTE1 Hepatocyte-specific data combined with Tao hepatocyte specific genes do not 
#        yield any matching genes (without secondary data.)
#-------------------------------------------------------------------------------------#

# filter for gene name
# See NOTE 1
# lapply(race.all, function(x) filter(x, gene_name %in% SCdata)) -> hepatocyte.all

#-------------------------------------------------------------------------------------#
# Section: Merge Single-Cell Hepatocyte data with race gene data object (2/2)
# + Search Single-cell hepatocyte-specific data secondary alias for gene name
#-------------------------------------------------------------------------------------#

# duplicate data frame where gene-name is duplicated for every alias name
SCdata$alias <- as.character(SCdata$alias)
# unnnest gene alias names that are separated by ','
SCdata %>% unnest(alias = strsplit(alias, ","), keep_empty = TRUE) -> unnestedSCdata

# drop gene_name column and use alias for matching genes, rename alias column to gene_name
unnestedSCdata %>% select(-gene_name) %>% rename(gene_name = alias) -> unnestedSCdata.dif

# filter for alias name for all races
lapply(race.all, function(x) inner_join(x, unnestedSCdata.dif, by  = "gene_name")) -> alias.filt

# retain and write to csv
names(alias.filt) <- races[1:3]
lapply(1:length(alias.filt), function(i) write.csv(alias.filt[[i]],
                                                  file = paste0("./02_output/06_SChepatocyteOutput/", names(alias.filt[i]), "_SChepatocyte_retained_stage1.csv"), row.names = FALSE))

# join all races for overlap
joinfunction(alias.filt) -> hepatocyte_allrace.joined

# write to csv
write.csv(hepatocyte_allrace.joined, paste0("./02_output/06_SChepatocyteOutput/allrace_SChepatocyte_retained_stage1.csv"))