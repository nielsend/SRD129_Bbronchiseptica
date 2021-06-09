#######################################################################
#SRD129 16S alpha and beta diversity
#Kathy Mou

#NOTES: 
#This code analyzes alpha and beta diversity statistics for nasal samples BB and Control groups only, 
#no days -12 and -6, and associated NMDS plots
#This script uses files created in "SRD129_phyloseq_each_tissue3.R"

#Clear workspace and load necessary packages
rm(list=ls())

#Set working directory on Mac desktop:
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
library(scales)
library(RColorBrewer)
library(philentropy)
library(cowplot)
library('wesanderson')
library("ggsci")

#Load functions at the very bottom of this script
NMDS_ellipse
veganCovEllipse

#Load image file
load("SRD129_alpha_beta_diversity3.RData")



#Save image file
save.image(file="SRD129_alpha_beta_diversity3.RData")

###########################################################################################################


########################################### Nasal no bad samples ######################################

nw.meta.nobad <- data.frame(phyloseq.nw.nobad@sam_data) #make phyloseq.nw.nobad sam_data into dataframe
nw.otu.nobad <- data.frame(t(phyloseq.nw.nobad@otu_table)) #make phyloseq.nw.nobad otu_table into dataframe
class(nw.meta.nobad) #data.frame
rownames(nw.meta.nobad) == row.names(nw.otu.nobad) #make sure rownames between nw.meta.nobad and nw.otu.nobad match exactly. Yes they do.
nw.meta.nobad$numOTUS <- rowSums(nw.otu.nobad > 1)
head(nw.meta.nobad)
dim(nw.meta.nobad) #460 6


#Sum of how many samples in each day*treatment group:
nrow(nw.meta.nobad[nw.meta.nobad$Day == "DNEG12",]) #40
nrow(nw.meta.nobad[nw.meta.nobad$Day == "DNEG6",]) #39
nrow(nw.meta.nobad[nw.meta.nobad$Day == "D0",]) #40
nrow(nw.meta.nobad[nw.meta.nobad$Day == "D1",]) #39
nrow(nw.meta.nobad[nw.meta.nobad$Day == "D3",]) #39
nrow(nw.meta.nobad[nw.meta.nobad$Day == "D7",]) #39
nrow(nw.meta.nobad[nw.meta.nobad$Day == "D10",]) #40
nrow(nw.meta.nobad[nw.meta.nobad$Day == "D14",]) #40
nrow(nw.meta.nobad[nw.meta.nobad$Day == "D21",]) #39
nrow(nw.meta.nobad[nw.meta.nobad$Day == "D28",]) #27
nrow(nw.meta.nobad[nw.meta.nobad$Day == "D36",]) #39
nrow(nw.meta.nobad[nw.meta.nobad$Day == "D42",]) #39


#Removing BB, PRRSV, DNEG12, DNEG6 from nw.meta
nw.meta2 <- nw.meta.nobad[!grepl("IAV", nw.meta.nobad$Treatment),]
nw.meta2 <- nw.meta2[!grepl("PRRSV", nw.meta2$Treatment),]
nw.meta2 <- nw.meta2[!grepl("DNEG12", nw.meta2$Day),]
nw.meta2 <- nw.meta2[!grepl("DNEG6", nw.meta2$Day),]
unique(nw.meta2$Day) #D0  D10 D14 D21 D28 D3  D36 D42 D7  D1 
unique(nw.meta2$Treatment) #BB      Control


#Sum of how many samples in each day*treatment group after removing bad samples, IAV, PRRSV groups
nrow(nw.meta2[nw.meta2$Day == "D0",]) #20
nrow(nw.meta2[nw.meta2$Day == "D1",]) #19
nrow(nw.meta2[nw.meta2$Day == "D3",]) #19
nrow(nw.meta2[nw.meta2$Day == "D7",]) #19
nrow(nw.meta2[nw.meta2$Day == "D10",]) #20
nrow(nw.meta2[nw.meta2$Day == "D14",]) #20
nrow(nw.meta2[nw.meta2$Day == "D21",]) #20
nrow(nw.meta2[nw.meta2$Day == "D28",]) #14
nrow(nw.meta2[nw.meta2$Day == "D36",]) #20
nrow(nw.meta2[nw.meta2$Day == "D42",]) #20

#Merge nw.meta2 with nw.otu and edit data frames
class(nw.otu.nobad) #dataframe
dim(nw.otu.nobad) #460 2196
dim(nw.otu) #472 2196
class(nw.meta2) #dataframe
dim(nw.meta2) #191 6
nw.otu2 <- merge(nw.meta2, nw.otu.nobad, by='row.names')
dim(nw.otu2) #191 2196
colnames(nw.otu2)
rownames(nw.otu2)
rownames(nw.otu2) <- nw.otu2$Row.names
nw.otu2[1] <- NULL
nrow(nw.otu2) #191
ncol(nw.otu2) #2196
head(nw.otu2[,1:6])
nw.otu2 <- nw.otu2[,-c(1:6)]
head(nw.otu2[2190:2196])
tail(nw.otu2, n=5)
nrow(nw.otu2) #191
ncol(nw.otu2) #2196
head(nw.meta2)
rownames(nw.meta2)
rownames(nw.otu2)
rownames(nw.meta2) == row.names(nw.otu2) #make sure rownames between two files match exactly. Yes they do.
nw.meta2$numOTUS <- rowSums(nw.otu2 > 1)
head(nw.meta2)
dim(nw.meta2) #191 9
colnames(nw.meta2)
rownames(nw.meta2)
rownames(nw.otu2)
unique(nw.meta2$Day) #D0  D10 D14 D21 D28 D3  D36 D42 D7  D1
unique(nw.meta2$Treatment) #BB Control

###NMDS calculation and plot k=2
nw.otu2[1:10,1:10]
##run NMDS_ellipse and veganCovEllipse functions
nw_NMDS.2 <- NMDS_ellipse(nw.meta2, nw.otu2, grouping_set = 'All')
#grouping_set: group by "All"
#run the NMDS_ellipse function (at end of script) to add to R functions
#NMDS_ellipse function created by Jules
#also run veganCovEllipse function as well

#Result: [1] "Ordination stress: 0.16599863901607"

#separating meta data and ellipse data to two lists to make NMDS plot
head(nw_NMDS.2)
nw_NMDS.2[[1]]
head(nw_NMDS.2[[1]])
tail(nw_NMDS.2[[1]])
nw.metanmds.2 <- nw_NMDS.2[[1]] 
#nw.metanmds.nobad has meta data + MDS calculations. Select this 1st list of nw_NMDS.2 using double brackets
nw.df_ell.2 <- nw_NMDS.2[[2]]
#nw.df_ell.nobad is accessing 2nd list from nw_NMDS.nobad that has ellipse calculations
#need two separate lists for making NMDS plot
nw.df_ell.2$group #20 levels
head(nw.df_ell.2)

#need Day and Treatment column for faceting
nw.df_ell.2$Day <- sub('_[A-Za-z]+', '\\1', nw.df_ell.2$group) #create Day column
#'\\1' Return the first part of the regular expression D[A-z0-9]+ from nw.df_ell$group

#Testing regular expressions
#thing1 <- sub('D[0-9]+', '\\1', 'D0_BB')
#thing2 <- sub('_[A-Za-z]+', '\\1', 'D0_BB')
#thing2

nw.df_ell.2$Treatment <- sub('D[A-z0-9]+_', '\\2', nw.df_ell.2$group) #create Treatment column
#'\\2' Return the second part of the sub expression ([A-Za-z]+) from nw.df_ell$group
head(nw.df_ell.2)
tail(nw.df_ell.2)
unique(nw.df_ell.2$Treatment) #BB Control

####NMDS plot for nasal
#Day only
(ggplot(data=nw.metanmds.2, aes(x=MDS1, y=MDS2, color=Day)) + geom_point() + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = .5) + 
  geom_path(data=nw.df_ell.2, aes(x=NMDS1, y=NMDS2, color=Day, group=group)) + 
  labs(caption = 'Ordination stress = 0.17') + scale_colour_brewer(palette="Paired"))

#Treatment only
ggplot(data=nw.metanmds.2, aes(x=MDS1, y=MDS2, color=Treatment)) + geom_point() + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = .5) + 
  geom_path(data=nw.df_ell.2, aes(x=NMDS1, y=NMDS2, color=Treatment, group=group)) + 
  #scale_color_brewer(palette="Dark2") +
  labs(caption = 'Ordination stress = 0.17')  

#Restructure level for nw.metanmds.2 and nw.df_ell.2
nw.metanmds.2$Day = factor(nw.metanmds.2$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
nw.df_ell.2$Day = factor(nw.df_ell.2$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
levels(nw.df_ell.2$Day)
levels(nw.metanmds.2$Day)

#Day and treatment (gray background)
ggplot(data=nw.metanmds.2, aes(x=MDS1, y=MDS2, color=Treatment)) + geom_point() + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = 0.5) + 
  geom_path(data=nw.df_ell.2, aes(x=NMDS1, y=NMDS2, color=Treatment, group=group)) + 
  facet_wrap(~Day, scales = 'free', nrow=2) +
  #scale_color_brewer(palette="Dark2") +
  theme_gray(base_size = 10) +
  theme(strip.text.x = element_text(size=15)) +
  labs(color="Treatment group")

#Create nw.metanmds.2a object (for gray points of the All days and treatments plot)
nw.metanmds.2a <- nw.metanmds.2
nw.metanmds.2a$Treatment2 = nw.metanmds.2a$Treatment
nw.metanmds.2a$Treatment2 <- as.character(nw.metanmds.2a$Treatment2)

#All days and treatments faceted by day (gridlines)
ggplot(nw.metanmds.2, aes(x=MDS1, y=MDS2)) +  annotate(x=nw.metanmds.2a$MDS1, y=nw.metanmds.2a$MDS2, color='grey57', geom = 'point')+
  geom_path(data = nw.df_ell.2, aes(x=NMDS1, y=NMDS2, color=Treatment), size=1.25) + 
  geom_point(aes(color = Treatment), size = 2) + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=Treatment), alpha=.5) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.border = element_rect(fill = NA, color = 'grey57'),
        axis.line = element_blank()) + facet_wrap(~Day, nrow = 2) +
  theme_bw() +
  labs(caption = 'Ordination stress = 0.17', color="Treatment group")
#Too difficult to see


#All days and treatments faceted by day (no gridlines or axes)
ggplot(nw.metanmds.2, aes(x=MDS1, y=MDS2)) +  annotate(x=nw.metanmds.2a$MDS1, y=nw.metanmds.2a$MDS2, color='grey57', geom = 'point')+
  geom_path(data = nw.df_ell.2, aes(x=NMDS1, y=NMDS2, color=Treatment), size=1.25) + 
  geom_point(aes(color = Treatment), size = 2) + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=Treatment), alpha=.5) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.border = element_rect(fill = NA, color = 'grey57'),
        axis.line = element_blank()) + facet_wrap(~Day, nrow = 2) +
  labs(subtitle = 'Ordination stress = 0.17')
#Too difficult to see




#Plotting just control group, with gridlines and axes, gray points for All days
#Subset control group
nw.metanmds.2acontrol <- subset(nw.metanmds.2a, Treatment== "Control", select = c(Pig:centroidY))
nw.metanmds.2control <- subset(nw.metanmds.2, Treatment == 'Control', select = c(Pig:centroidY))
nw.df_ell.2control <- subset(nw.df_ell.2, Treatment=='Control', select =c(NMDS1:Treatment))

#Plot control group, gridlines, axes, gray points for All days
ggplot(nw.metanmds.2control, aes(x=MDS1, y=MDS2)) +  annotate(x=nw.metanmds.2acontrol$MDS1, y=nw.metanmds.2acontrol$MDS2, color='grey57', geom = 'point')+
  geom_path(data = nw.df_ell.2control, aes(x=NMDS1, y=NMDS2, color=Treatment), size=1.25) + 
  geom_point(aes(color = Treatment), size = 2) + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=Treatment), alpha=.5) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.border = element_rect(fill = NA, color = 'grey57'),
        axis.line = element_blank()) + facet_wrap(~Day, nrow = 2) +
  theme_bw() +
  labs(caption = 'Ordination stress = 0.16')











######################---------BETA DIVERSITY STATS AND PLOTS WITH PAIRWISE ADONIS: NASAL---------######################

#Nasal no bad samples
nw.adon.2 <- pairwise.adonis(nw.otu2, nw.meta2$All, sim.method = 'bray', p.adjust.m = 'bonferroni')
nw.adon.2
write.csv(nw.adon.2, file='SRD129_Nasal_AllPairwiseAdonisBBControl.csv', row.names=TRUE)
nw.adon.2$pairs #list all comparisons in the "pairs" column
goodcomps.nw.2 <- c(grep('D0_[A-Za-z]+ vs D0_[A-Za-z]+', nw.adon.2$pairs),
                  grep('D1_[A-Za-z]+ vs D1_[A-Za-z]+', nw.adon.2$pairs),
                  grep('D3_[A-Za-z]+ vs D3_[A-Za-z]+', nw.adon.2$pairs),
                  grep('D7_[A-Za-z]+ vs D7_[A-Za-z]+', nw.adon.2$pairs),
                  grep('D10_[A-Za-z]+ vs D10_[A-Za-z]+', nw.adon.2$pairs),
                  grep('D14_[A-Za-z]+ vs D14_[A-Za-z]+', nw.adon.2$pairs),
                  grep('D21_[A-Za-z]+ vs D21_[A-Za-z]+', nw.adon.2$pairs),
                  grep('D28_[A-Za-z]+ vs D28_[A-Za-z]+', nw.adon.2$pairs),
                  grep('D36_[A-Za-z]+ vs D36_[A-Za-z]+', nw.adon.2$pairs),
                  grep('D42_[A-Za-z]+ vs D42_[A-Za-z]+', nw.adon.2$pairs))
#use regular expressions
#[A-Za-z] match all capital, lowercase letters
#+ match a whole word and not just one letter (if you didn't have "+")
#c creates the vector, lump all pairs of specific interest groups together
#make vector to combine all the pairwise comparisons you're interested in (same day, different treatment group)
goodcomps.nw.2
nw.adon.good.2 <- nw.adon.2[goodcomps.nw.2,] #rename this vector
nw.adon.good.2
nw.adon.good.2$p.adjusted <- p.adjust(nw.adon.good.2$p.value, method = 'fdr')
nw.adon.good.2$p.adjusted2 <- round(nw.adon.good.2$p.adjusted, 3)
nw.adon.good.2$p.adjusted2[nw.adon.good.2$p.adjusted2 > 0.05] <- NA
nw.adon.good.2
write.csv(nw.adon.good.2, file='SRD129_Nasal_SelectPairwiseAdonisBBControlNoDNEG.csv', row.names=TRUE)



########################---------ALPHA DIVERSITY STATS AND PLOTS WITH PAIRWISE WILCOX: NASAL---------###################

#######################################-----Nasal no bad samples------######################################################

####average shannon, invsimpson, numOTUs per "All" subtype for nobad samples
nw.meta2$shannon <- diversity(nw.otu2) #diversity is vegan function, default index is shannon
#added a shannon index column in nw.meta.nobad
nw.meta2$invsimpson <- diversity(nw.otu2,index = 'invsimpson')
#use invsimpson since easier to understand than simpson
#don't need to "inverse" simpson values
#simpson: low the number = higher diversity

levels(sample_data(nw.meta2)$Day) #"D0"  "D1"  "D10" "D14" "D21" "D28" "D3"  "D36" "D42" "D7"
nw.meta2$Day = factor(nw.meta2$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
levels(sample_data(nw.meta2)$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"   
head(nw.meta2)

####average shannon, invsimpson, numOTUs per "All" subtype
nw.2.average.shannon.invsimpson.numOTUs <- aggregate(nw.meta2[, 6:8], list(nw.meta2$All), mean)
print(nw.2.average.shannon.invsimpson.numOTUs)

#         Group.1   numOTUS  shannon invsimpson
#1        D0_BB  44.50000 2.057835   4.566800
#2   D0_Control  71.10000 2.553987   6.469838
#3        D1_BB  48.33333 2.270197   5.944310
#4   D1_Control  50.70000 2.265464   5.130755
#5       D10_BB  73.10000 2.671551   7.719212
#6  D10_Control  61.50000 2.441225   5.103775
#7       D14_BB  74.30000 2.502470   4.178867
#8  D14_Control  80.20000 2.821534   6.861513
#9       D21_BB  85.90000 3.089974   8.074825
#10 D21_Control  78.70000 2.938500   7.259222
#11      D28_BB  97.60000 3.309728  11.343234
#12 D28_Control 102.75000 3.460608  13.326738
#13       D3_BB  86.44444 2.942788   7.719057
#14  D3_Control  82.90000 2.910306   7.380954
#15      D36_BB  92.20000 3.408973  14.059936
#16 D36_Control  90.10000 3.146692  13.230399
#17      D42_BB 131.70000 3.923845  29.600788
#18 D42_Control 111.00000 3.521438  17.944710
#19       D7_BB  65.66667 2.464188   5.793326
#20  D7_Control  56.80000 2.378387   5.480581

####box and whisker plot: shannon
#Renaming injected and infeed to parenteral and medicated feed
#nw.meta$Treatment2 = nw.meta$Treatment
#nw.meta$Treatment2 <- as.character(nw2.meta$Treatment2)
#nw2.meta$Treatment2[nw2.meta$Treatment2 == 'Injected'] <- "Parenteral"
#nw2.meta$Treatment2[nw2.meta$Treatment2 == 'Infeed'] <- "Medicated Feed"
nw.shan.bw.2 <- ggplot(data = nw.meta2, aes(x=Treatment, y=shannon, group=All, fill=Treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~Day, scales = 'free') +
  scale_y_continuous(name = "Shannon diversity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text.x = element_text(size=10),
        axis.title.y = element_text(size=15))
#free allows each plot to customize scale to the specific data set (no forced scale applied to all plots)
#position = position_dodge2(preserve = 'total') --> fixes the ggplot box width;
#prevents narrow boxes from forming in plot (they become wider)
nw.shan.bw.2
#To change treatment group colors
#nw.shan <- nw.shan + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "none")


####box and whisker plot: inverse simpson
nw.invsimp.bw.2 <- ggplot(data = nw.meta2, aes(x=Treatment, y=invsimpson, group=All, fill=Treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~Day, scales = 'free') +
  scale_y_continuous(name = "Inverse Simpson diversity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text.x = element_text(size=10),
        axis.title.y = element_text(size=15))
nw.invsimp.bw.2
#To change treatment group colors
#nw.invsimp.bw <- nw.invsimp.bw + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "none")

#combining plots
nw.combalpha.2 <- plot_grid(nw.shan.bw.2 + theme(legend.position = "none"), nw.invsimp.bw.2, labels = "AUTO")
nw.combalpha.2


####pairwise wilcox test shannon

#day 0
meta0.nw.2 <- nw.meta2[nw.meta2$Day == 'D0',] 
pairwise.wilcox.test(meta0.nw.2$shannon, meta0.nw.2$All, p.adjust.method = 'none') 
#Pairwise comparisons using Wilcoxon rank sum test 
#data:  meta0.nw.2$shannon and meta0.nw.2$All 
#           D0_BB
#D0_Control 0.043

#day 1
meta1.nw.2 <- nw.meta2[nw.meta2$Day == 'D1',] 
pairwise.wilcox.test(meta1.nw.2$shannon, meta1.nw.2$All, p.adjust.method = 'none') 
#Pairwise comparisons using Wilcoxon rank sum test 
#data:  meta1.nw.2$shannon and meta1.nw.2$All 

#           D1_BB
#D1_Control 0.5    

#day 3
meta3.nw.2 <- nw.meta2[nw.meta2$Day == 'D3',]
pairwise.wilcox.test(meta3.nw.2$shannon, meta3.nw.2$All, p.adjust.method = 'none')
#              D3_BB
#   D3_Control 0.84

#day 7
meta7.nw.2 <- nw.meta2[nw.meta2$Day == 'D7',]
pairwise.wilcox.test(meta7.nw.2$shannon, meta7.nw.2$All, p.adjust.method = 'none')
#                D7_BB
#     D7_Control 0.6 

#day 10
meta10.nw.2 <- nw.meta2[nw.meta2$Day == 'D10',]
pairwise.wilcox.test(meta10.nw.2$shannon, meta10.nw.2$All, p.adjust.method = 'none')
#                   D10_BB
#       D10_Control 0.31  

#day 14
meta14.nw.2 <- nw.meta2[nw.meta2$Day == 'D14',]
pairwise.wilcox.test(meta14.nw.2$shannon, meta14.nw.2$All, p.adjust.method = 'none')
#                 D14_BB
#     D14_Control 0.25 

#day 21
meta21.nw.2 <- nw.meta2[nw.meta2$Day == 'D21',]
pairwise.wilcox.test(meta21.nw.2$shannon, meta21.nw.2$All, p.adjust.method = 'none')
#                 D21_BB
#     D21_Control 0.22

#day 28
meta28.nw.2 <- nw.meta2[nw.meta2$Day == 'D28',]
pairwise.wilcox.test(meta28.nw.2$shannon, meta28.nw.2$All, p.adjust.method = 'none')
#                 D28_BB
#     D28_Control 0.54 

#day 36
meta36.nw.2 <- nw.meta2[nw.meta2$Day == 'D36',]
pairwise.wilcox.test(meta36.nw.2$shannon, meta36.nw.2$All, p.adjust.method = 'none')
#                 D36_BB
#     D36_Control 0.53

#day 42
meta42.nw.2 <- nw.meta2[nw.meta2$Day == 'D42',]
pairwise.wilcox.test(meta42.nw.2$shannon, meta42.nw.2$All, p.adjust.method = 'none')
#                 D42_BB
#     D42_Control 0.17 

typeof(nw.meta2) #list
dim(nw.meta2) #191 8
head(nw.meta2[,1:8]) #shannon column at the end
nw.2.pairwise.wilcox.shannon.test <- pairwise.wilcox.test(nw.meta2$shannon, nw.meta2$All, p.adjust.method = 'none')
print(nw.2.pairwise.wilcox.shannon.test)
options(max.print=1000000)
#get results of all pairwise wilcox test comparisons
#copied printout to excel






####pairwise wilcox test inverse simpson

#day 0
meta0.nw.inv.2 <- nw.meta2[nw.meta2$Day == 'D0',]
pairwise.wilcox.test(meta0.nw.inv.2$invsimpson, meta0.nw.inv.2$All, p.adjust.method = 'none')
#                D0_BB
#     D0_Control 0.17

#day 1
meta1.nw.inv.2 <- nw.meta2[nw.meta2$Day == 'D1',]
pairwise.wilcox.test(meta1.nw.inv.2$invsimpson, meta1.nw.inv.2$All, p.adjust.method = 'none')
#                   D1_BB
#       D1_Control  0.97 

#day 3
meta3.nw.inv.2 <- nw.meta2[nw.meta2$Day == 'D3',]
pairwise.wilcox.test(meta3.nw.inv.2$invsimpson, meta3.nw.inv.2$All, p.adjust.method = 'none')
#                   D3_BB
#       D3_Control  0.84

#day 7
meta7.nw.inv.2 <- nw.meta2[nw.meta2$Day == 'D7',]
pairwise.wilcox.test(meta7.nw.inv.2$invsimpson, meta7.nw.inv.2$All, p.adjust.method = 'none')
#            D7_BB
#     D7_Control 0.84

#day 10
meta10.nw.inv.2 <- nw.meta2[nw.meta2$Day == 'D10',]
pairwise.wilcox.test(meta10.nw.inv.2$invsimpson, meta10.nw.inv.2$All, p.adjust.method = 'none')
#                 D10_BB
#     D10_Control 0.22 

#day 14
meta14.nw.inv.2 <- nw.meta2[nw.meta2$Day == 'D14',]
pairwise.wilcox.test(meta14.nw.inv.2$invsimpson, meta14.nw.inv.2$All, p.adjust.method = 'none')
#                 D14_BB
#     D14_Control 0.052 

#day 21
meta21.nw.inv.2 <- nw.meta2[nw.meta2$Day == 'D21',]
pairwise.wilcox.test(meta21.nw.inv.2$invsimpson, meta21.nw.inv.2$All, p.adjust.method = 'none')
#                 D21_BB
#     D21_Control 0.39

#day 28
meta28.nw.inv.2 <- nw.meta2[nw.meta2$Day == 'D28',]
pairwise.wilcox.test(meta28.nw.inv.2$invsimpson, meta28.nw.inv.2$All, p.adjust.method = 'none')
#                 D28_BB
#     D28_Control 0.64

#day 36
meta36.nw.inv.2 <- nw.meta2[nw.meta2$Day == 'D36',]
pairwise.wilcox.test(meta36.nw.inv.2$invsimpson, meta36.nw.inv.2$All, p.adjust.method = 'none')
#                 D36_BB
#     D36_Control 0.48     

#day 42
meta42.nw.inv.2 <- nw.meta2[nw.meta2$Day == 'D42',]
pairwise.wilcox.test(meta42.nw.inv.2$invsimpson, meta42.nw.inv.2$All, p.adjust.method = 'none')
#                 D42_BB
#     D42_Control 0.19   

nw.2.pairwise.wilcox.invsimpson.test <- pairwise.wilcox.test(nw.meta2$invsimpson, nw.meta2$All, p.adjust.method = 'none')
#get results of all pairwise wilcox test comparisons for inverse simpson
print(nw.2.pairwise.wilcox.invsimpson.test)









######################################--------------FUNCTIONS---------------#####################################
#NMDS functions
#When I changed the distance_method to "jaccard", I still got the same ordination stress
#and same plots
NMDS_ellipse <- function(metadata, OTU_table, grouping_set, distance_method = 'bray', rand_seed = 77777, MDS_trymax = 1000){
  generic_dist <- vegdist(OTU_table, method = distance_method)
  set.seed(rand_seed)
  generic_MDS <- metaMDS(generic_dist, k = 2, trymax = MDS_trymax, autotransform = FALSE)
  stress <- generic_MDS$stress
  nmds_points <- as.data.frame(generic_MDS$points)
  
  metadata <- metadata[match(rownames(generic_MDS$points), rownames(metadata)),]
  
  metanmds <- cbind(metadata, nmds_points)
  
  ord <- ordiellipse(generic_MDS, metanmds[[grouping_set]], label = TRUE, conf = .95, kind = 'se', draw = 'none')
  
  nmds.mean <- aggregate(metanmds[,grep("MDS", colnames(metanmds))], list(group=metanmds[[grouping_set]]), mean)
  
  metanmds[[grouping_set]] <- factor(metanmds[[grouping_set]]) # this 'set' needs to be passed in from function
  
  #check to make sure at least 3 obs for each grouping_set
  
  numobs <- metanmds %>% group_by_(grouping_set) %>% summarise(n=n())
  if(!(all(is.element(c(0,1,2), numobs$n)))){
    
    df_ell <- data.frame()
    for (d in levels(metanmds[[grouping_set]])){
      df_ell <- rbind(df_ell, cbind(as.data.frame(with(metanmds[metanmds[[grouping_set]] == d,],
                                                       veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
    }
    
    
    levels(metanmds[[grouping_set]])
    
    # this loop assigns the group centroid X coordinates to each sample
    metanmds$centroidX <- NA
    metanmds$centroidY <- NA
    
    
    for (level in levels(metanmds[[grouping_set]])){
      metanmds[metanmds[[grouping_set]] == level,]$centroidX <- nmds.mean$MDS1[nmds.mean$group == level]
      metanmds[metanmds[[grouping_set]] == level,]$centroidY <- nmds.mean$MDS2[nmds.mean$group == level]
      
      
    }
    print(paste('Ordination stress:', stress, sep = ' '))
    return(list(metanmds, df_ell))
  }
  print('One of your groups in "grouping_set" has less than 3 observations, cannot generate elipses')
  
  
}

veganCovEllipse <- function (cov, center = c(0,0), scale = 1, npoints = 100){
  # this is a function in the vegan package?
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}






#Copied from researchgate: https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}
