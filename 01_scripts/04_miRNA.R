#-------------------------------------------------------------------------------------#
# Project: Biomarker identification
# Purpose: Identify biomarkers from literature in BG output for Stage A
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