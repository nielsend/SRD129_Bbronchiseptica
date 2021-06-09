#######################################################################
#SRD129 Nasal Tonsil Magnitude of Change
#Kathy Mou

#NOTES: 
#This code plots the magnitude of change for nasal and tonsil beta diversity

#Clear workspace and load necessary packages
rm(list=ls())

#Set working directory (either on desktop or network drive)
#Mac
setwd("~/Desktop/Microbiome/Projects/SRD129/SRD129_2000singletons")
setwd("~/Desktop/SRD129/SRD129_2000singletons")


#Load library packages
library(vegan)
library(ggplot2)
#ggplot:  fill = color in boxes
#         color = colored outline
#aes: aesthetics declared
#alt - to get <- symbol
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(philentropy)
library(cowplot)
library(ggrepel)

#Load image file
load("SRD129_NasalTonsil_MagnitudeOfChange3.RData")





#Save image file
save.image(file="SRD129_NasalTonsil_MagnitudeOfChange3.RData")

######################################################################

#Nasal
nasalmag <- read.csv("SRD129_Nasal_BBMagnitudeofChange.csv") #used F values from SRD129_Nasal_SelectPairwiseAdonis_nobad.xlsx
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

#BB and control only
nasalmag3<- nasalmag[!grepl("PRRSV", nasalmag$Treatment),]
nasalmag3<- nasalmag3[!grepl("IAV", nasalmag3$Treatment),]
nasalmag3<- nasalmag3[!grepl("-12", nasalmag3$Day),]
nasalmag3<- nasalmag3[!grepl("-6", nasalmag3$Day),]
unique(nasalmag3$Day) #0  1  3  7  10 14 21 28 36 42
unique(nasalmag3$Treatment) #BB
nasalmag3
(nasalmag3a <- ggplot(data=nasalmag3, aes(x=Day, y=Value, group=Treatment)) +
  geom_line(aes(color=Treatment)) + geom_point(aes(color=Treatment)) +
  ylab("PERMANOVA F vs control \n(difference relative to control)") +
  scale_fill_manual(values=c("blue")) +
  scale_color_manual(values=c("blue")) +
  geom_label_repel(aes(Day, Value, label=Value), box.padding = 0.35, point.padding=0.5, segment.color = 'grey50'))


#Too much text
#nasalmag2 <- nasalmag2 + geom_label_repel(aes(Day, Value, label=Value), box.padding = 0.35, point.padding=0.5,
#      segment.color = 'grey50') +
#      theme_classic(base_size = 12)
#nasalmag2

#Tonsil
tonsilmag <- read.csv("SRD129_TonsilMagnitudeofChange.csv")
class(tonsilmag) #data.frame
tonsilmag$Day <- factor(tonsilmag$Day)
tonsilmag$Value <- round(tonsilmag$Value, 2)
tonsilmag2 <- ggplot(data=tonsilmag, aes(x=Day, y=Value, group=Treatment)) +
  geom_line(aes(color=Treatment)) + geom_point(aes(color=Treatment)) +
  ylab("PERMANOVA F vs control \n(difference relative to control)") +
  scale_fill_manual(values=c("green", "blue", "purple")) +
  scale_color_manual(values=c("green", "blue", "purple"))
tonsilmag2

#Too much text
#tonsilmag2 <- tonsilmag2 + geom_label_repel(aes(Day, Value, label=Value), box.padding = 0.35, point.padding=0.5,
#                             segment.color = 'grey50') +
#  theme_classic(base_size = 12)
#tonsilmag2

#Combining plots
nwtt.combmag <- plot_grid(nasalmag2, tonsilmag2, labels = "AUTO")
nwtt.combmag
