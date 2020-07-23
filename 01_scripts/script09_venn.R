#-------------------------------------------------------------------------------------#
# Project: Venn Diagram
# Purpose: Visualize number of genes from each race for Stage 1
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
library(VennDiagram)

#-------------------------------------------------------------------------------------#
# Input files
# + Load BG output for transcripts
#-------------------------------------------------------------------------------------#
rm(list = ls())
inputpaths <- dir("./00_data/stage1", full.names = TRUE)
race_dfs <- lapply(inputpaths, function(x) read.delim(x))
races <- c("asian", "black", "white")
names(race_dfs) <- races

#-------------------------------------------------------------------------------------#
# Functions
# + clean differentially expressed genes
# + retain top FC differentially expressed genes
# + merge race df's function
# + differential expression fold change cutoff
#-------------------------------------------------------------------------------------#

clean_deg <- function(df, cutoff){
  df %>% select(gene_name, fc, qval, pval) %>%
    mutate(logfc = log2(fc)) %>% mutate(abs_fc = abs(logfc)) %>% 
    filter(qval < 0.05 & abs_fc > cutoff) %>% 
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
  ABWmerge %>% 
    rename(pval.white = pval) %>%
    rename(abs_fc.white = abs_fc) %>%
    rename(fc.white = fc) %>%
    rename(logfc.white = logfc) %>%
    rename(qval.white =  qval) -> ABWmerge
  print(paste0("# genes for Asian: ", nrow(listobject[[1]])))
  print(paste0("# genes for Black: ", nrow(listobject[[2]])))
  print(paste0("# genes for White: ", nrow(listobject[[3]])))
  print(paste0("Total number of overlapped genes: ", nrow(ABWmerge)))
  return(ABWmerge)
}

#-------------------------------------------------------------------------------------#
# Section: Clean data
# + create necessary data columns with 'clean_deg'
# + filter for unique transcripts by largest FC for each transcript
# + filter for absolute fold change > 1.0
#-------------------------------------------------------------------------------------#

race_dfs <- lapply(race_dfs, clean_deg, 1.0)

race_dfs <- lapply(race_dfs, retain_topfc)



#-------------------------------------------------------------------------------------#
# Create venn diagram
# + overlap DEGs for each race
# + set venn variables to overlapped DEGs
#-------------------------------------------------------------------------------------#

area1 = nrow(race_dfs$asian)
area2 = nrow(race_dfs$black)
area3 = nrow(race_dfs$white)
n12 = nrow(inner_join(race_dfs$asian, race_dfs$black, by = "gene_name"))
n23 = nrow(inner_join(race_dfs$black, race_dfs$white, by = "gene_name"))
n13 = nrow(inner_join(race_dfs$asian, race_dfs$white, by = "gene_name"))
n123 = nrow(joinfunction(race_dfs))

# create grid for venn diagram
grid.newpage()

# draw venn diagram 
draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12   = n12, n23 = n23, n13 = n13,
                 n123 = n123, 
                 category = c(paste0("Asian (", area1, ")"),
                            paste0("African-American (", area2,")"),
                            paste0("Caucasian (", area3, ")")),
                 lty      = "blank", 
                 fill     = c("blue3", "pink3", "green4"),
                 alpha = c(0.5, 0.5, 0.3), 
                 cat.fontface = "bold",
                 fontfamily = "Arial",
                 cat.fontfamily = "Arial",
                 cat.pos = c(335, 15, 180),
                 cat.just = list(c(-0.5, 0),  # asian
                                 c(0.5,1),    # white
                                 c(0.5,0.3)), # black
                 cat.cex = 1.5,  # size for each category name
                 cex = 2.5,  # size for each area label
                 units = px,
                 height = 3000,
                 width = 3000,
                 resolution = 800,
                 output = TRUE,
                 margin = 0.025)

