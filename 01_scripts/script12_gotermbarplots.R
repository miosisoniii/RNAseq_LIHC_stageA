#-------------------------------------------------------------------------------------#
# Project: TCGA-LIHC
# Purpose: Barplot Visualizations for GOterms
# Author: Artemio Sison III
# R Version: 4.0.1 "See Things Now"
#-------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------#
# Load Dependencies
#-------------------------------------------------------------------------------------#
rm(list = ls())
library(dplyr)
library(ggplot2)
library(cowplot)


#-------------------------------------------------------------------------------------#
# Read GO term tables
#-------------------------------------------------------------------------------------#
# read in all names
inputpaths <- dir("../10_goterm", full.names = TRUE)
race_dfs <- lapply(inputpaths, function(x) read.csv(x))

# create names for each dataframe in race_dfs
races <- c("asian", "black", "white")

# assign names to list of race dataframes
names(race_dfs) <- races[1:3]

# view dataframe
asian <- race_dfs$asian
#-------------------------------------------------------------------------------------#
# Clean data
#-------------------------------------------------------------------------------------#

# create function for cleaning data
cleanGO <- function(data) {
  data %>% 
    select(X..background.genes, X..genes, category, 
           description, FDR.value, term.name) %>%
    rename(bg.genes = X..background.genes,
           genes = X..genes) %>%
    filter(category == "GO Component" | category == "GO Process" | category == "GO Function")
}


# use function and add race name
asian <- cleanGO(race_dfs$asian) %>% mutate(race = "asian")
black <- cleanGO(race_dfs$black) %>% mutate(race = "black")
white <- cleanGO(race_dfs$white) %>% mutate(race = "white")

# use reduce() 
# create list of df's
combinedGOlist <- list(asian, black, white)

# nested merge
combinedGO <- merge(merge(
  combinedGOlist[[1L]],
  combinedGOlist[[2L]], all = TRUE),
  combinedGOlist[[3L]], all = TRUE)

# filter for rows with GO terms that contain 3 occurrences (all 3 races)
combinedGO <- combinedGO %>% group_by(description) %>% filter(n() == 3)

# group by race 
topGO <- combinedGO %>% arrange(FDR.value, .by_group = TRUE)

#-------------------------#
# Calculate mean FDR
#-------------------------#
# calculate mean FDR
aveGOterm <- aggregate(FDR.value ~ description, data = topGO, FUN = mean)

# identify top FDR 
ave_GOterm <- top_n(aveGOterm, 25)

# match names of top FDR within topGO
topGOout <- topGO %>% filter(description %in% ave_GOterm$description)

# recode FDR to calculate log FDR
topGOout <- topGOout %>% mutate(FDR.value = -log(FDR.value))

# Add colors in conditional statement
topGOtest <- topGOout %>% mutate(color = if_else(category == "GO Function", "green",
                                                 if_else(category == "GO Component", "blue", "red")))




#-------------------------------------------------------------------------------------#
# Testing Barplot Code
#-------------------------------------------------------------------------------------#
# This code arranges the order of the race color in each bar
topGOtesting <- topGOtest[order(topGOtest$race),]


topGOtesting1 <- topGOtesting
# create specific order to match the geom text
topGOtesting1$race <- with(topGOtesting1, reorder(race, -FDR.value))


plot <- ggplot(topGOtesting, aes(fill = race, 
                              group = category,
                              # Arranges bars in ascending order
                              y = reorder(description, -FDR.value),
                              x = FDR.value)) 

plot + 
  theme(axis.text = element_text(size = 12),
        axis.text.y = element_text(color = topGOtesting$color)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(aes(label = color, col = color)) 
#   theme_set(theme_classic()+
#               theme(panel.grid.major.x = element_line(colour = 'grey60', linetype = 'dashed'),
#                     panel.grid.major.y = element_line(colour = 'grey60', linetype = 'dashed'),
#                     axis.ticks.y = element_blank(),
#                     axis.text.x = element_text(colour = 'black', size = 16),
#                     axis.ticks.x = element_line(colour = 'grey60'),
#                     axis.ticks.length = unit(3, "mm"),
#                     aspect.ratio = (600/450),
#                     axis.title.x=element_blank(),
#                     axis.title.y=element_blank()))
#   
# 
# topGOtesting1 <- topGOtesting %>% arrange(-FDR.value) %>% 
#   mutate(colour = factor(color, color), 
#          description = factor(description, description))
# 
# 
# 
# 
























# reorder data for plotting
# topGOtest <- topGOtest[order(topGOtest$race, decreasing = FALSE),] # original code
topGOtest <- topGOtest[order(topGOtest$race),]
 
#-------------------------------------------------------------------------------------#
# Plot Stacked Barplot
#-------------------------------------------------------------------------------------#
# https://stackoverflow.com/questions/50768148/ggplot-retaining-axis-label-coloring-with-reordered-data
topGOtest$description <- with(topGOtest, reorder(description, -FDR.value))
cols <- c("GO Process" = "red", "GO Function" = "green", "GO Component" = "blue")

# a <- ifelse(ave_GOtermtestcolor$category == "GO Function", "green",
#             ifelse(ave_GOtermtestcolor$category == "GO Component", "blue", "red"))
# a <- ifelse(topGOtest$category == "GO Function", "green",
#             ifelse(topGOtest$category == "GO Component", "blue", "red"))

plot <- ggplot(topGOtest, aes(fill = race, 
                              group = category,
                              y = description,
                              #y = reorder(description, -FDR.value),
                              x = FDR.value))

plot + 
  # cowplot theme
  #theme_minimal_vgrid() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_bar(position = "stack", stat = "identity") +
  # original code that worked
  # theme(axis.text = element_text(size = 12),
  #        axis.text.y = element_text(color = ave_GOtermtestcolor$color)) +
  theme(axis.text = element_text(size = 12),
         axis.text.y = element_text(color = topGOtest$color)) +
        # legend.text = element_text(color = cols[topGOtest$color])) +
  
  # theme(axis.text.y = element_text(color = ifelse(topGOout$category == "GO Process", "red",
  #                                                 ifelse(topGOout$category == "GO Function", "green", "blue")))) +
  # theme(axis.text.y = element_text(color = ifelse(ave_GOtermtestcolor$category == "GO Process", "red",
  #                                                 ifelse(ave_GOtermtestcolor$category == "GO Function", "green", "blue")))) +

  # legend
  # scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"),
  #                   name = "GO Term\nCategory",
  #                   # breaks = c("GO Process", "GO Function", "GO Component"),
  #                   breaks = c("asian", "black", "white"),
  #                   labels = c("Asian", "African-\nAmerican", "Caucasian")) +
scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"),
                  name = "Race",
                  # breaks = c("GO Process", "GO Function", "GO Component"),
                  breaks = c("asian", "black", "white"),
                  labels = c("Asian", "African-\nAmerican", "Caucasian")) +
  
  # labels
  labs(title = "Gene Ontology Term Representation Across Race and GO Classification", 
       y = "Gene Ontology Term",
       x = "-Log FDR Value",
       fill = "Race") +
  
  # this can be used for reference
  # http://zevross.com/blog/2014/08/04/beautiful-plotting-in-r-a-ggplot2-cheatsheet-3/#manually-adding-legend-items-guides-override.aes
  
  # This adds a label to the plot, but I can use this to map a legend
  # remove "a" from geom legend
  # https://stackoverflow.com/questions/18337653/remove-a-from-legend-when-using-aesthetics-and-geom-text
  geom_text(aes(label = category, col = category)) +
  scale_color_manual(name="GO Term\nCategory", values = c("GO Process" = "red", "GO Function" = "green", "GO Component" = "blue"))

# Add code for positioning of legend to inside of plot


  
  
ggsave("../02_output/11_barplot/stackedbarplot_draft.jpg", width = 12, height = 8, units = "in", dpi = 800)
 






#-------------------------#
# Create mapping for based on average
#-------------------------#
# 
# # select only relevant columns - desdcription (key), category
# sel_GOterm <- topGOtest %>% select(description, category) %>% distinct(description, .keep_all = TRUE)
# 
# # add category based on description
# ave_GOtermtest <- inner_join(ave_GOterm, sel_GOterm, by = "description")
# 
# # color based on catgoery
# ave_GOtermtestcolor <- ave_GOtermtest %>% mutate(color = if_else(category == "GO Function", "green",
#                                                                  if_else(category == "GO Component", "blue", "red")))
# 
# # order using same order from previous code for entire dataset on mapping df
# ave_GOtermtestcolor$description <- with(ave_GOtermtestcolor, reorder(description, log(FDR.value)))