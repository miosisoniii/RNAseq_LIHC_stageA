#-------------------------------------------------------------------------------------#
# Project: TCGA LIHC Paper
# Purpose: Create Barplots for Overlapping Genes
# Author: Artemio Sison III
# R Version: 4.0.1 "See Things Now"
#-------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------#
# Load Dependencies
#-------------------------------------------------------------------------------------#
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
rm(list = ls())

#-------------------------------------------------------------------------------------#
# Load input data
#-------------------------------------------------------------------------------------#
data <- read.csv("../02_output/02_overlapDEGoutput/stage1/allrace_overlapDEG_stage1_FC1.csv")

#-------------------------------------------------------------------------------------#
# Clean data
#-------------------------------------------------------------------------------------#
dataselect <- data %>% select(gene_name, logfc.asian, logfc.black, logfc.white)

dataselect <- dataselect %>% gather(race, "logfc", 2:4)

# use dplyr::separate to separate the columns!
dataselect <- dataselect %>% 
  separate(col = race, into = c("drop", "race"), sep = "(\\.)") %>%
  select(-drop)
#-------------------------------------------------------------------------------------#
# Plot
#-------------------------------------------------------------------------------------#
plot <- geom_bar(data = dataselect, aes(x = gene_name, y = ))