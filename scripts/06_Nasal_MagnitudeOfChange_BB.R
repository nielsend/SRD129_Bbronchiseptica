#######################################################################
#SRD129 16S - Magnitude of Change in Nasal Microbiota
#By Mou, KT

#NOTES: 
#This code plots the magnitude of change in nasal microbiota of BB group relative to Control group

#Clear workspace and load necessary packages
rm(list=ls())

#File needed:
#bb.pairwisecomparisons.csv: copy the F.Model, R2, p.value, p.adjusted, p.adjusted2 columns to a new spreadsheet as columns #3-7. 
#For the first column, label as "Day" and list all days (0, 1, 3, 7, 10, etc.). Label second column as "Treatment" and list "BB" in all respective cells under this column. 
#Save this as SRD129_BB_MagnitudeOfChange.csv.

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

######################################################################

#Nasal
nasalmag <- read.csv("./data/SRD129_BB_MagnitudeOfChange.csv") #created 
class(nasalmag) #data.frame
nasalmag$Day <- factor(nasalmag$Day)
nasalmag
nasal <- ggplot(data=nasalmag, aes(x=Day, y=F.Model, group=Treatment)) +
  geom_line(aes(color=Treatment)) + geom_point(aes(color=Treatment)) +
  ylab("PERMANOVA F vs control \n(difference relative to control)") +
  scale_fill_manual(values=c("#F8766D")) +
  scale_color_manual(values=c("#F8766D")) +  
  geom_label_repel(aes(label=p.adjusted2), box.padding = 0.35, point.padding=0.5,segment.color = 'grey50') +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(size = 14), 
        axis.text.x = element_text(size=14), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14), 
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14)) +
  labs(color="Treatment group")
nasal

#Save 'nasal' as a .tiff for publication, 300dpi
ggsave("Figure_2.tiff", plot=nasal, width = 8, height = 4, dpi = 300, units =c("in"))