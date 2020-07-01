#-------------------------------------------------------------------------------------#
# Project: TCGA LIHC Stage A
# Purpose: Use API to perform differential expression of Child Pugh Stage A samples
# Author: Artemio Sison III
# R Version: 4.0.1 "See Things Now"
# Date Created: 20200615
# Notes: https://www.bioconductor.org/packages/release/bioc/vignettes/sevenbridges/inst/doc/bioc-workflow.html
# Authentication Code: XXXXX
#-------------------------------------------------------------------------------------#

# Dependencies
library("sevenbridges")

# Set environment for API using Authentication Token
sbg_set_env("https://cgc-api.sbgenomics.com/v2", "XXXX")

# Create Auth object using credentials from environment
a <- Auth(from = "env")

#-------------------------------------------------------------------------------------#
# Section: Interact with API to download files
# Description: This code downloads the transcript output files for viewing
# ################################################################################### #
# + Assign project to variable "p"
# + Filter for BG output files (ie. ballgown transcript list)
# + Download files
#-------------------------------------------------------------------------------------#


# List files belonging to project
p <- a$project(id = "rachelzayas_pi/tcga-lihc-disparities")
#head(p$file())

# Filter for Asian/White/Black transcript FPKM files
# Found in http://docs.cancergenomicscloud.org/docs/get-your-authentication-token
# Section 3.9.2.4

# Filter files: 
p$file(
  #metadata = list(
    # race = "asian",
    # race = "black or african american",
    # race = "white"
  #),
  tag = c("bg_stage1out")
) -> transcripts

# Move filtered files into working R Environment via downloading
# Download transcripts to ./00_data folder
download(transcripts, "./00_data/stage1")
