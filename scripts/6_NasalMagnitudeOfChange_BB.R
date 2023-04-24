#######################################################################
#SRD129 16S - Magnitude of Change in Nasal Microbiota


#Purpose: Plot the F-statistic from PERMANOVA pairwise comparisons of Control and BB groups, over time (displays the magnitude of change in the nasal bacterial community structure of the BB group relative to control)

#Files needed:
#BB_MagnitudeOfChange.csv

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(philentropy)
library(cowplot)
library(ggrepel)


#######################################################################

#The F-values were obtained from 'BB.pairwisecomparisons.csv' data generated in the "Beta diversity" section
#BB_MagnitudeOfChange.csv was created from BB.pairwisecomparisons.csv
#To make the BB_MagnitudeOfChange.csv file, open BB.pairwisecomparisons.csv file in excel, 
#copy columns "F. Model" through "p.adjusted2" and paste in a separate spreadsheet. Add "Day" and "Treatment" columns and save as "BB_MagnitudeOfChange.csv".
View(BB.pairwisecomparisons)

BB.mag <- read.csv("./BB_MagnitudeOfChange.csv")
colnames(BB.mag) <- c("F.Model","R2","p.value","p.adjusted","p.adjusted2","Day","Treatment")
class(BB.mag)
BB.mag$Day <- factor(BB.mag$Day) #Encode "Day" as a factor
BB.mag$Day = factor(BB.mag$Day, levels=c("D0", "D1", "D3","D7", "D10","D14", "D21", "D36", "D42"))

BB.mag2 <- ggplot(data=BB.mag, aes(x=Day, y=F.Model, group=Treatment)) +
  geom_line(aes(color=Treatment), size=0.7) + geom_point(aes(color=Treatment)) +
  ylab("PERMANOVA F vs control \n(difference relative to control)") +
  scale_fill_manual(values=c("#F8766D")) +
  scale_color_manual(values=c("#F8766D")) +
  geom_label_repel(aes(label=p.adjusted), box.padding = 0.35, point.padding=1,segment.color = 'grey50') +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size=14), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), legend.text=element_text(size=14), legend.title=element_text(size=14)) +
  labs(color="Treatment group")
BB.mag2

#Save 'BB.mag2' as a .tiff for publication, 500dpi
ggsave("MagOfChange.tiff", plot=BB.mag2, width = 10, height = 5, dpi = 500, units =c("in"))
