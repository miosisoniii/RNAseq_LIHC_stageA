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
plot <- ggplot(data = dataselect, aes(fill = race, y = logfc, x = gene_name))

# Stacked Barplot
# plot + geom_bar(position = "stack", stat = "identity")
  
# Dodged barplot
plot + geom_bar(position = "dodge", stat = "identity")



#-------------------------------------------------------------------------------------#
#  Things to be added to the plots
#-------------------------------------------------------------------------------------#

# Change gene name to transcripts
# Labels on the bottom
# Change logFC
# Change title to Stage I
# Labels 

# colors
colors <- c("blue", "orange", "gray")


# Dodged barplot
final <- plot + geom_bar(position = "dodge", stat = "identity") +
  # Labels for LogFC on each bar, adjusted for each bar label
  geom_text(aes(label = round(logfc, 2)), position = position_dodge(width = 0.9), vjust = 1.6) +
  ggtitle("Stage A HCC Liver Tissue Overlapping Transcripts") +
  xlab("Transcript") +
  ylab(expression(Log[2]~"Fold Change"))

# Save final image with appropriate size/dimensions
ggsave("../02_output/11_barplot/stageAoverlap.jpeg", width = 8, height = 8)
