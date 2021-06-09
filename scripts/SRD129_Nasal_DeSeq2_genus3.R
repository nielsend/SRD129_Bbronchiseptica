#####################################################################################################
#SRD129 Nasal DESeq2 - Genus
#Kathy Mou

#NOTES: 
#This code analyzes differential abundance of OTUs in nasal samples between BB and control groups
#using DESeq2 - GENUS ONLY, no bad samples, no DNEG12, DNEG6
#also log2fold change and basemean plots

#Clear workspace and load necessary packages
rm(list=ls())

#Set working directory (either on local or network drive, Mac or PC)
#Mac
setwd("~/Desktop/Microbiome/Projects/SRD129/SRD129_2000singletons")
setwd("C:/Users/Kathy.Mou/Desktop/SRD129/SRD129_2000singletons")


#Load library packages

library(DESeq2)
library(phyloseq)
library(ggplot2)
library(tidyr)
library("wesanderson")
library(plotly)
library(gapminder)
library("ggsci")

#Load image file
load("SRD129_Nasal_DESeq2_genus3.RData")




#Save image file
save.image(file="SRD129_Nasal_DESeq2_genus3.RData")

#Additional notes
#Jules says plots describes the log-fold changes seen in differential abundance plots as
#enriched for "x" taxa in "x group"... but why negative log-fold change for one group and
#positive fold change for the other group?
#You can add genus to your plots, and use those as supplemental figures

#Annotations
#BC = Bordetella, control
#IC = IAV, control
#PC = PRRSV, control

#Set plots to have gray background with white gridlines
theme_set(theme_gray())

########################################################################################################

####### PREPARING OBJECTS FOR DESEQ2 ANALYSIS ########

#Load files
otu2 <- import_mothur(mothur_shared_file = 'SRD129.outsingletons.abund.opti_mcc.shared') #use unrarified data
taxo2 <- import_mothur(mothur_constaxonomy_file = 'SRD129.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
taxo2
meta2 <- read.table(file = 'SRD129metadata06292018nobad.csv', sep = ',', header = TRUE) #no bad samples

#Organize meta file
rownames(meta2) <- meta2$Sample
meta2 <- meta2[,-2] #remove Sample column
meta2 <- meta2[,-1] #remove X column
meta2 <- meta2[,-5] #remove All column and replace with new Set column (below)
#meta2$Day <- gsub("D", "", meta2$Day) # remove "D"
meta2$Set <- paste(meta2$Day, meta2$Tissue, meta2$Treatment, sep = '_')

#Make phyloseq object SRD129 (combine taxonomy, OTU, and metadata)
phy_meta2 <- sample_data(meta2) 
SRD129 <- phyloseq(otu2, taxo2)
SRD129 <- merge_phyloseq(SRD129, phy_meta2)   # combines the metadata with this phyloseq object
colnames(tax_table(SRD129))
colnames(tax_table(SRD129)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
SRD129

#Prune
SRD129 <- prune_samples(sample_sums(SRD129) > 2000, SRD129)  # This removes samples that have fewer than 2000 sequences associated with them.
#Jules used 1000. You use 2000 as a way to subsample.
SRD129 <- prune_taxa(taxa_sums(SRD129) > 10, SRD129)        # removes OTUs that occur less than 10 times globally
tax_table(SRD129) [1:5, 1:6] #see what's in tax_table first 5 rows, first 6 columns

# If you want to group OTUs uncomment the tax_glom() line and select your desired taxrank
# right now all these analysis are done at the OTU level.

SRD129.genus <- tax_glom(SRD129, taxrank = "Genus")
# This method merges species that have the same taxonomy at a certain taxanomic rank. 
# Its approach is analogous to tip_glom, but uses categorical data instead of a tree. 

#Subset nasal samples to count number of sequences for each OTU
SRD129.nw <- subset_samples(SRD129.genus, Tissue == 'N')
SRD129.nw.1 <- psmelt(SRD129.nw)
write.csv(SRD129.nw.1, file="SRD129_GenusSampleData.csv")

SRD129.nw.otu <- otu_table(SRD129.nw)
write.csv(SRD129.nw.otu, file="SRD129_DESeqOTUtable.csv")

######################################## PRIMARY COMPARISONS TO MAKE ####################################################

# NMDS plot showed that disperion is different between days, so subset by day and tissue

# NMDS plot of important comparisons to make between treatment and control (significant changes in beta diversity between treatments)

# Day 0, 1, 3: PRRSV
# Day 7, 10: BB, PRRSV
# Day -6, 14, 21, 28, 36, 42: BB, PRRSV, IAV

# NMDS plot of other comparisons to make (no significant changes in beta diversity between treatments)



##################################################################################################################################

###############################################BB comparison between days###############################################
sample_data(SRD129.genus)
SRD129.BB.nw <- subset_samples(SRD129.genus, Treatment == 'BB' & Tissue == 'N')
sample_sums(SRD129.BB.nw)
colnames(otu_table(SRD129.BB.nw)) #check on all the sample names
SRD129.BB.nw <- prune_taxa(taxa_sums(SRD129.BB.nw) > 1, SRD129.BB.nw)
#if taxa_sums is >1, then it will print that out in SRD129.DNEG12.nw object and not include anything with <1.
rowSums(SRD129.BB.nw@otu_table)
SRD129.BB.nw.De <- phyloseq_to_deseq2(SRD129.BB.nw, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
SRD129.BB.nw.De <- DESeq(SRD129.BB.nw.De, test = "Wald", fitType = "parametric")


###############################################Comparing Day 0 to day 1###############################################
meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D0_N_BB")
#BB = 10
sum(meta2$Set == "D1_N_BB")
#BB = 9

#Extract results from a DESeq analysis, organize table
res.D0_1.bb = results(SRD129.BB.nw.De, contrast=c("Set","D1_N_BB","D0_N_BB"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D0_1.bb
sigtab.D0_1.bb = res.D0_1.bb[which(res.D0_1.bb$padj < .05), ]
sigtab.D0_1.bb = cbind(as(sigtab.D0_1.bb, "data.frame"), as(tax_table(SRD129.BB.nw)[rownames(sigtab.D0_1.bb), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D0_1.bb$padj, scientific = TRUE)
sigtab.D0_1.bb$newp <- format(round(sigtab.D0_1.bb$padj, digits = 3), scientific = TRUE)
sigtab.D0_1.bb$Treatment <- ifelse(sigtab.D0_1.bb$log2FoldChange >=0, "Day_1", "Day_0")
#ifelse: Assigning "Day_1" = yes, "Day_0" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D0_1.bb

#Summarize sigtab.D0.bc
sum.sigtab.D0_1.bb <- summary(sigtab.D0_1.bb)
sum.sigtab.D0_1.bb

#ggplot
deseq.D0_1.bb <- ggplot(sigtab.D0_1.bb, aes(x=reorder(rownames(sigtab.D0_1.bb), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0_1.bb), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group on Day 1 Relative to Day 0 at Nasal Site')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Day_1", "Day_0"), values = c('#E69F00', '#999999'))
deseq.D0_1.bb

#Add OTU and comparisons columns
sigtab.D0_1.bb$OTU <- rownames(sigtab.D0_1.bb)
sigtab.D0_1.bb$comp <- 'BB_nasal_D1vsD0'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab.bb <- sigtab.D0_1.bb



###############################################Comparing Day 0 to day 3###############################################
meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D0_N_BB")
#BB = 10
sum(meta2$Set == "D3_N_BB")
#BB = 9

#Extract results from a DESeq analysis, organize table
res.D0_3.bb = results(SRD129.BB.nw.De, contrast=c("Set","D3_N_BB","D0_N_BB"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D0_3.bb
sigtab.D0_3.bb = res.D0_3.bb[which(res.D0_3.bb$padj < .05), ]
sigtab.D0_3.bb = cbind(as(sigtab.D0_3.bb, "data.frame"), as(tax_table(SRD129.BB.nw)[rownames(sigtab.D0_3.bb), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D0_3.bb$padj, scientific = TRUE)
sigtab.D0_3.bb$newp <- format(round(sigtab.D0_3.bb$padj, digits = 3), scientific = TRUE)
sigtab.D0_3.bb$Treatment <- ifelse(sigtab.D0_3.bb$log2FoldChange >=0, "Day_3", "Day_0")
#ifelse: Assigning "Day_3" = yes, "Day_0" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D0_3.bb

#Summarize sigtab.D0.bc
sum.sigtab.D0_3.bb <- summary(sigtab.D0_3.bb)
sum.sigtab.D0_3.bb

#ggplot
deseq.D0_3.bb <- ggplot(sigtab.D0_3.bb, aes(x=reorder(rownames(sigtab.D0_3.bb), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0_3.bb), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group on Day 3 Relative to Day 0 at Nasal Site')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Day_3", "Day_0"), values = c('#E69F00', '#999999'))
deseq.D0_3.bb

#Add OTU and comparisons columns
sigtab.D0_3.bb$OTU <- rownames(sigtab.D0_3.bb)
sigtab.D0_3.bb$comp <- 'BB_nasal_D3vsD0'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab.bb <- rbind(final.nonsigtab.bb,sigtab.D0_3.bb)



###############################################Comparing Day 0 to day 7###############################################
meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D0_N_BB")
#BB = 10
sum(meta2$Set == "D7_N_BB")
#BB = 9

#Extract results from a DESeq analysis, organize table
res.D0_7.bb = results(SRD129.BB.nw.De, contrast=c("Set","D7_N_BB","D0_N_BB"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D0_7.bb
sigtab.D0_7.bb = res.D0_7.bb[which(res.D0_7.bb$padj < .05), ]
sigtab.D0_7.bb = cbind(as(sigtab.D0_7.bb, "data.frame"), as(tax_table(SRD129.BB.nw)[rownames(sigtab.D0_7.bb), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D0_7.bb$padj, scientific = TRUE)
sigtab.D0_7.bb$newp <- format(round(sigtab.D0_7.bb$padj, digits = 3), scientific = TRUE)
sigtab.D0_7.bb$Treatment <- ifelse(sigtab.D0_7.bb$log2FoldChange >=0, "Day_7", "Day_0")
#ifelse: Assigning "Day_7" = yes, "Day_0" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D0_7.bb

#Summarize sigtab.D0.bc
sum.sigtab.D0_7.bb <- summary(sigtab.D0_7.bb)
sum.sigtab.D0_7.bb

#ggplot
deseq.D0_7.bb <- ggplot(sigtab.D0_7.bb, aes(x=reorder(rownames(sigtab.D0_7.bb), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0_7.bb), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group on Day 7 Relative to Day 0 at Nasal Site')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Day_7", "Day_0"), values = c('#E69F00', '#999999'))
deseq.D0_7.bb

#Add OTU and comparisons columns
sigtab.D0_7.bb$OTU <- rownames(sigtab.D0_7.bb)
sigtab.D0_7.bb$comp <- 'BB_nasal_D7vsD0'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab.bb <- rbind(final.nonsigtab.bb,sigtab.D0_7.bb)




###############################################Comparing Day 0 to day 10###############################################
meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D0_N_BB")
#BB = 10
sum(meta2$Set == "D10_N_BB")
#BB = 10

#Extract results from a DESeq analysis, organize table
res.D0_10.bb = results(SRD129.BB.nw.De, contrast=c("Set","D10_N_BB","D0_N_BB"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D0_10.bb
sigtab.D0_10.bb = res.D0_10.bb[which(res.D0_10.bb$padj < .05), ]
sigtab.D0_10.bb = cbind(as(sigtab.D0_10.bb, "data.frame"), as(tax_table(SRD129.BB.nw)[rownames(sigtab.D0_10.bb), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D0_10.bb$padj, scientific = TRUE)
sigtab.D0_10.bb$newp <- format(round(sigtab.D0_10.bb$padj, digits = 3), scientific = TRUE)
sigtab.D0_10.bb$Treatment <- ifelse(sigtab.D0_10.bb$log2FoldChange >=0, "Day_10", "Day_0")
#ifelse: Assigning "Day_10" = yes, "Day_0" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D0_10.bb

#Summarize sigtab.D0.bc
sum.sigtab.D0_10.bb <- summary(sigtab.D0_10.bb)
sum.sigtab.D0_10.bb

#ggplot
deseq.D0_10.bb <- ggplot(sigtab.D0_10.bb, aes(x=reorder(rownames(sigtab.D0_10.bb), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0_10.bb), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group on Day 10 Relative to Day 0 at Nasal Site')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Day_10", "Day_0"), values = c('#E69F00', '#999999'))
deseq.D0_10.bb

#Add OTU and comparisons columns
sigtab.D0_10.bb$OTU <- rownames(sigtab.D0_10.bb)
sigtab.D0_10.bb$comp <- 'BB_nasal_D10vsD0'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab.bb <- rbind(final.nonsigtab.bb,sigtab.D0_10.bb)



###############################################Comparing Day 0 to day 14###############################################
meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D0_N_BB")
#BB = 10
sum(meta2$Set == "D14_N_BB")
#BB = 10

#Extract results from a DESeq analysis, organize table
res.D0_14.bb = results(SRD129.BB.nw.De, contrast=c("Set","D14_N_BB","D0_N_BB"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D0_14.bb
sigtab.D0_14.bb = res.D0_14.bb[which(res.D0_14.bb$padj < .05), ]
sigtab.D0_14.bb = cbind(as(sigtab.D0_14.bb, "data.frame"), as(tax_table(SRD129.BB.nw)[rownames(sigtab.D0_14.bb), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D0_14.bb$padj, scientific = TRUE)
sigtab.D0_14.bb$newp <- format(round(sigtab.D0_14.bb$padj, digits = 3), scientific = TRUE)
sigtab.D0_14.bb$Treatment <- ifelse(sigtab.D0_14.bb$log2FoldChange >=0, "Day_14", "Day_0")
#ifelse: Assigning "Day_14" = yes, "Day_0" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D0_14.bb

#Summarize sigtab.D0.bc
sum.sigtab.D0_14.bb <- summary(sigtab.D0_14.bb)
sum.sigtab.D0_14.bb

#ggplot
deseq.D0_14.bb <- ggplot(sigtab.D0_14.bb, aes(x=reorder(rownames(sigtab.D0_14.bb), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0_14.bb), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group on Day 14 Relative to Day 0 at Nasal Site')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Day_14", "Day_0"), values = c('#E69F00', '#999999'))
deseq.D0_14.bb

#Add OTU and comparisons columns
sigtab.D0_14.bb$OTU <- rownames(sigtab.D0_14.bb)
sigtab.D0_14.bb$comp <- 'BB_nasal_D14vsD0'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab.bb <- rbind(final.nonsigtab.bb,sigtab.D0_14.bb)



###############################################Comparing Day 0 to day 21###############################################
meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D0_N_BB")
#BB = 10
sum(meta2$Set == "D21_N_BB")
#BB = 10

#Extract results from a DESeq analysis, organize table
res.D0_21.bb = results(SRD129.BB.nw.De, contrast=c("Set","D21_N_BB","D0_N_BB"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D0_21.bb
sigtab.D0_21.bb = res.D0_21.bb[which(res.D0_21.bb$padj < .05), ]
sigtab.D0_21.bb = cbind(as(sigtab.D0_21.bb, "data.frame"), as(tax_table(SRD129.BB.nw)[rownames(sigtab.D0_21.bb), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D0_21.bb$padj, scientific = TRUE)
sigtab.D0_21.bb$newp <- format(round(sigtab.D0_21.bb$padj, digits = 3), scientific = TRUE)
sigtab.D0_21.bb$Treatment <- ifelse(sigtab.D0_21.bb$log2FoldChange >=0, "Day_21", "Day_0")
#ifelse: Assigning "Day_21" = yes, "Day_0" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D0_21.bb

#Summarize sigtab.D0.bc
sum.sigtab.D0_21.bb <- summary(sigtab.D0_21.bb)
sum.sigtab.D0_21.bb

#ggplot
deseq.D0_21.bb <- ggplot(sigtab.D0_21.bb, aes(x=reorder(rownames(sigtab.D0_21.bb), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0_21.bb), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group on Day 21 Relative to Day 0 at Nasal Site')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Day_21", "Day_0"), values = c('#E69F00', '#999999'))
deseq.D0_21.bb

#Add OTU and comparisons columns
sigtab.D0_21.bb$OTU <- rownames(sigtab.D0_21.bb)
sigtab.D0_21.bb$comp <- 'BB_nasal_D21vsD0'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab.bb <- rbind(final.nonsigtab.bb,sigtab.D0_21.bb)


###############################################Comparing Day 0 to day 28###############################################
meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D0_N_BB")
#BB = 10
sum(meta2$Set == "D28_N_BB")
#BB = 10

#Extract results from a DESeq analysis, organize table
res.D0_28.bb = results(SRD129.BB.nw.De, contrast=c("Set","D28_N_BB","D0_N_BB"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D0_28.bb
sigtab.D0_28.bb = res.D0_28.bb[which(res.D0_28.bb$padj < .05), ]
sigtab.D0_28.bb = cbind(as(sigtab.D0_28.bb, "data.frame"), as(tax_table(SRD129.BB.nw)[rownames(sigtab.D0_28.bb), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D0_28.bb$padj, scientific = TRUE)
sigtab.D0_28.bb$newp <- format(round(sigtab.D0_28.bb$padj, digits = 3), scientific = TRUE)
sigtab.D0_28.bb$Treatment <- ifelse(sigtab.D0_28.bb$log2FoldChange >=0, "Day_28", "Day_0")
#ifelse: Assigning "Day_28" = yes, "Day_0" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D0_28.bb

#Summarize sigtab.D0.bc
sum.sigtab.D0_28.bb <- summary(sigtab.D0_28.bb)
sum.sigtab.D0_28.bb

#ggplot
deseq.D0_28.bb <- ggplot(sigtab.D0_28.bb, aes(x=reorder(rownames(sigtab.D0_28.bb), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0_28.bb), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group on Day 28 Relative to Day 0 at Nasal Site')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Day_28", "Day_0"), values = c('#E69F00', '#999999'))
deseq.D0_28.bb

#Add OTU and comparisons columns
sigtab.D0_28.bb$OTU <- rownames(sigtab.D0_28.bb)
sigtab.D0_28.bb$comp <- 'BB_nasal_D28vsD0'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab.bb <- rbind(final.nonsigtab.bb,sigtab.D0_28.bb)



###############################################Comparing Day 0 to day 36###############################################
meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D0_N_BB")
#BB = 10
sum(meta2$Set == "D36_N_BB")
#BB = 10

#Extract results from a DESeq analysis, organize table
res.D0_36.bb = results(SRD129.BB.nw.De, contrast=c("Set","D36_N_BB","D0_N_BB"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D0_36.bb
sigtab.D0_36.bb = res.D0_36.bb[which(res.D0_36.bb$padj < .05), ]
sigtab.D0_36.bb = cbind(as(sigtab.D0_36.bb, "data.frame"), as(tax_table(SRD129.BB.nw)[rownames(sigtab.D0_36.bb), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D0_36.bb$padj, scientific = TRUE)
sigtab.D0_36.bb$newp <- format(round(sigtab.D0_36.bb$padj, digits = 3), scientific = TRUE)
sigtab.D0_36.bb$Treatment <- ifelse(sigtab.D0_36.bb$log2FoldChange >=0, "Day_36", "Day_0")
#ifelse: Assigning "Day_36" = yes, "Day_0" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D0_36.bb

#Summarize sigtab.D0.bc
sum.sigtab.D0_36.bb <- summary(sigtab.D0_36.bb)
sum.sigtab.D0_36.bb

#ggplot
deseq.D0_36.bb <- ggplot(sigtab.D0_36.bb, aes(x=reorder(rownames(sigtab.D0_36.bb), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0_36.bb), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group on Day 36 Relative to Day 0 at Nasal Site')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Day_36", "Day_0"), values = c('#E69F00', '#999999'))
deseq.D0_36.bb

#Add OTU and comparisons columns
sigtab.D0_36.bb$OTU <- rownames(sigtab.D0_36.bb)
sigtab.D0_36.bb$comp <- 'BB_nasal_D36vsD0'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab.bb <- rbind(final.nonsigtab.bb,sigtab.D0_36.bb)




###############################################Comparing Day 0 to day 42###############################################
meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D0_N_BB")
#BB = 10
sum(meta2$Set == "D42_N_BB")
#BB = 10

#Extract results from a DESeq analysis, organize table
res.D0_42.bb = results(SRD129.BB.nw.De, contrast=c("Set","D42_N_BB","D0_N_BB"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D0_42.bb
sigtab.D0_42.bb = res.D0_42.bb[which(res.D0_42.bb$padj < .05), ]
sigtab.D0_42.bb = cbind(as(sigtab.D0_42.bb, "data.frame"), as(tax_table(SRD129.BB.nw)[rownames(sigtab.D0_42.bb), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D0_42.bb$padj, scientific = TRUE)
sigtab.D0_42.bb$newp <- format(round(sigtab.D0_42.bb$padj, digits = 3), scientific = TRUE)
sigtab.D0_42.bb$Treatment <- ifelse(sigtab.D0_42.bb$log2FoldChange >=0, "Day_42", "Day_0")
#ifelse: Assigning "Day_42" = yes, "Day_0" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D0_42.bb

#Summarize sigtab.D0.bc
sum.sigtab.D0_42.bb <- summary(sigtab.D0_42.bb)
sum.sigtab.D0_42.bb

#ggplot
deseq.D0_42.bb <- ggplot(sigtab.D0_42.bb, aes(x=reorder(rownames(sigtab.D0_42.bb), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0_42.bb), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group on Day 42 Relative to Day 0 at Nasal Site')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Day_42", "Day_0"), values = c('#E69F00', '#999999'))
deseq.D0_42.bb

#Add OTU and comparisons columns
sigtab.D0_42.bb$OTU <- rownames(sigtab.D0_42.bb)
sigtab.D0_42.bb$comp <- 'BB_nasal_D42vsD0'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab.bb <- rbind(final.nonsigtab.bb,sigtab.D0_42.bb)

########################################################################################################################

write.csv(final.nonsigtab.bb, file= "SRD129_FinalDiffAbundNasalGenus_SignificantDaysWithinBBGroup.csv")

########################################################################################################################





########################################################################################################################

############################################################ Day -12 Nasal #######################################################

sample_data(SRD129.genus)

#Making DESeq object from phyloseq
SRD129.DNEG12.nw <- subset_samples(SRD129.genus, Day == 'DNEG12' & Tissue == 'N')
sample_sums(SRD129.DNEG12.nw)
colnames(otu_table(SRD129.DNEG12.nw)) #check on all the sample names
SRD129.DNEG12.nw <- prune_taxa(taxa_sums(SRD129.DNEG12.nw) > 1, SRD129.DNEG12.nw)
#if taxa_sums is >1, then it will print that out in SRD129.DNEG12.nw object and not include anything with <1.
rowSums(SRD129.DNEG12.nw@otu_table)
SRD129.DNEG12.nw.De <- phyloseq_to_deseq2(SRD129.DNEG12.nw, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
SRD129.DNEG12.nw.De <- DESeq(SRD129.DNEG12.nw.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing



######################################################### Day 0 Nasal #########################################################

sample_data(SRD129.genus)

SRD129.D0.nw <- subset_samples(SRD129.genus, Day == 'D0' & Tissue == 'N')
sample_sums(SRD129.D0.nw)
colnames(otu_table(SRD129.D0.nw)) #check on all the sample names
SRD129.D0.nw <- prune_taxa(taxa_sums(SRD129.D0.nw) > 1, SRD129.D0.nw)
#if taxa_sums is >1, then it will print that out in SRD129.D0.nw object and not include anything with <1.
rowSums(SRD129.D0.nw@otu_table)
SRD129.D0.nw.De <- phyloseq_to_deseq2(SRD129.D0.nw, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
SRD129.D0.nw.De <- DESeq(SRD129.D0.nw.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing


######### Day 0 Nasal BB vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D0_N_BB")
#BB = 10
sum(meta2$Set == "D0_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D0.bc = results(SRD129.D0.nw.De, contrast=c("Set","D0_N_BB","D0_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D0.bc
sigtab.D0.bc = res.D0.bc[which(res.D0.bc$padj < .05), ]
sigtab.D0.bc = cbind(as(sigtab.D0.bc, "data.frame"), as(tax_table(SRD129.D0.nw)[rownames(sigtab.D0.bc), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D0.bc$padj, scientific = TRUE)
sigtab.D0.bc$newp <- format(round(sigtab.D0.bc$padj, digits = 3), scientific = TRUE)
sigtab.D0.bc$Treatment <- ifelse(sigtab.D0.bc$log2FoldChange >=0, "BB", "Control")
#ifelse: Assigning "In-feed" = yes, "Control" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D0.bc

#Summarize sigtab.D0.bc
sum.sigtab.D0.bc <- summary(sigtab.D0.bc)
sum.sigtab.D0.bc

#ggplot
deseq.D0.bc <- ggplot(sigtab.D0.bc, aes(x=reorder(rownames(sigtab.D0.bc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.bc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group Relative to Control at Nasal Site on Day 0')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("BB", "Control"), values = c('#E69F00', '#999999'))
deseq.D0.bc

#Add OTU and comparisons columns
sigtab.D0.bc$OTU <- rownames(sigtab.D0.bc)
sigtab.D0.bc$comp <- 'D0_nasal_BBvsControl'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D0.bc)


######### Day 0 Nasal IAV vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D0_N_IAV")
#IAV = 10
sum(meta2$Set == "D0_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D0.nw.De$Set
res.D0.ic = results(SRD129.D0.nw.De, contrast=c("Set","D0_N_IAV","D0_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(SRD129.D0.nw.De, contrast=c("Set","D0_N_IAV","D0_N_Control")) 
sigtab.D0.ic = res.D0.ic[which(res.D0.ic$padj < .05), ]
sigtab.D0.ic = cbind(as(sigtab.D0.ic, "data.frame"), as(tax_table(SRD129.D0.nw)[rownames(sigtab.D0.ic), ], "matrix"))
format(sigtab.D0.ic$padj, scientific = TRUE)
sigtab.D0.ic$newp <- format(round(sigtab.D0.ic$padj, digits = 3), scientific = TRUE)
sigtab.D0.ic$Treatment <- ifelse(sigtab.D0.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D0.ic
sum.sigtab.D0.ic <- summary(sigtab.D0.ic)
sum.sigtab.D0.ic

#ggplot
deseq.D0.ic <- ggplot(sigtab.D0.ic, aes(x=reorder(rownames(sigtab.D0.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 0')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "IAV"), values = c('#999999', '#CC0066'))
deseq.D0.ic

#Add OTU and comparisons columns
sigtab.D0.ic
sigtab.D0.ic$OTU <- rownames(sigtab.D0.ic)
sigtab.D0.ic
sigtab.D0.ic$comp <- 'D0_nasal_IAVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D0.ic)


######### Day 0 Nasal PRRSV vs Control ###################

meta2$Set
#Number of pigs per group: 
sum(meta2$Set == "D0_N_PRRSV")
#PRRSV = 10
sum(meta2$Set == "D0_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D0.pc = results(SRD129.D0.nw.De, contrast=c("Set","D0_N_PRRSV","D0_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D0.pc
sigtab.D0.pc = res.D0.pc[which(res.D0.pc$padj < .05), ]
sigtab.D0.pc = cbind(as(sigtab.D0.pc, "data.frame"), as(tax_table(SRD129.D0.nw)[rownames(sigtab.D0.pc), ], "matrix"))
format(sigtab.D0.pc$padj, scientific = TRUE)
sigtab.D0.pc$newp <- format(round(sigtab.D0.pc$padj, digits = 3), scientific = TRUE)
sigtab.D0.pc$Treatment <- ifelse(sigtab.D0.pc$log2FoldChange >=0, "PRRSV", "Control")

#Summarize sigtab.D0.pc
sum.sigtab.D0.pc <- summary(sigtab.D0.pc)
sum.sigtab.D0.pc

#ggplot
deseq.D0.pc <- ggplot(sigtab.D0.pc, aes(x=reorder(rownames(sigtab.D0.pc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.pc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in PRRS Virus Group Relative to Control at Nasal Site on Day 0')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "PRRSV"), values = c('#999999', '#99FF66'))
deseq.D0.pc

#Add OTU and comparisons columns
sigtab.D0.pc
sigtab.D0.pc$OTU <- rownames(sigtab.D0.pc)
sigtab.D0.pc
sigtab.D0.pc$comp <- 'D0_nasal_PRRSVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D0.pc)







##################################################### Day 1 ######################################################################

############## Day 1 Nasal #########################

sample_data(SRD129.genus)
SRD129.D1.nw <- subset_samples(SRD129.genus, Day == 'D1' & Tissue == 'N')
sample_sums(SRD129.D1.nw)
colnames(otu_table(SRD129.D1.nw)) #check on all the sample names
SRD129.D1.nw <- prune_taxa(taxa_sums(SRD129.D1.nw) > 1, SRD129.D1.nw)
#if taxa_sums is >1, then it will print that out in SRD129.D1.nw object and not include anything with <1.
rowSums(SRD129.D1.nw@otu_table)
SRD129.D1.nw.De <- phyloseq_to_deseq2(SRD129.D1.nw, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
SRD129.D1.nw.De <- DESeq(SRD129.D1.nw.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing


######### Day 1 Nasal BB vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D1_N_BB")
#BB = 9
sum(meta2$Set == "D1_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D1.bc = results(SRD129.D1.nw.De, contrast=c("Set","D1_N_BB","D1_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D1.bc
sigtab.D1.bc = res.D1.bc[which(res.D1.bc$padj < .05), ]
sigtab.D1.bc = cbind(as(sigtab.D1.bc, "data.frame"), as(tax_table(SRD129.D1.nw)[rownames(sigtab.D1.bc), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D1.bc$padj, scientific = TRUE)
sigtab.D1.bc$newp <- format(round(sigtab.D1.bc$padj, digits = 3), scientific = TRUE)
sigtab.D1.bc$Treatment <- ifelse(sigtab.D1.bc$log2FoldChange >=0, "BB", "Control")
#ifelse: Assigning "In-feed" = yes, "Control" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D1.bc

#Summarize sigtab.D1.bc
sum.sigtab.D1.bc <- summary(sigtab.D1.bc)
sum.sigtab.D1.bc

#ggplot
deseq.D1.bc <- ggplot(sigtab.D1.bc, aes(x=reorder(rownames(sigtab.D1.bc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D1.bc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group Relative to Control at Nasal Site on Day 1')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("BB", "Control"), values = c('#E69F00', '#999999'))
deseq.D1.bc

#Add OTU and comparisons columns
sigtab.D1.bc$OTU <- rownames(sigtab.D1.bc)
sigtab.D1.bc$comp <- 'D1_nasal_BBvsControl'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D1.bc)


######### Day 1 Nasal IAV vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D1_N_IAV")
#IAV = 10
sum(meta2$Set == "D1_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D1.nw.De$Set
res.D1.ic = results(SRD129.D1.nw.De, contrast=c("Set","D1_N_IAV","D1_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(SRD129.D1.nw.De, contrast=c("Set","D1_N_IAV","D1_N_Control")) 
sigtab.D1.ic = res.D1.ic[which(res.D1.ic$padj < .05), ]
sigtab.D1.ic = cbind(as(sigtab.D1.ic, "data.frame"), as(tax_table(SRD129.D1.nw)[rownames(sigtab.D1.ic), ], "matrix"))
format(sigtab.D1.ic$padj, scientific = TRUE)
sigtab.D1.ic$newp <- format(round(sigtab.D1.ic$padj, digits = 3), scientific = TRUE)
sigtab.D1.ic$Treatment <- ifelse(sigtab.D1.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D1.ic
sum.sigtab.D1.ic <- summary(sigtab.D1.ic)
sum.sigtab.D1.ic

#ggplot
deseq.D1.ic <- ggplot(sigtab.D1.ic, aes(x=reorder(rownames(sigtab.D1.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D1.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 1')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "IAV"), values = c('#999999', '#CC0066'))
deseq.D1.ic

#Add OTU and comparisons columns
sigtab.D1.ic
sigtab.D1.ic$OTU <- rownames(sigtab.D1.ic)
sigtab.D1.ic
sigtab.D1.ic$comp <- 'D1_nasal_IAVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D1.ic)


######### Day 1 Nasal PRRSV vs Control ###################

meta2$Set
#Number of pigs per group: 
sum(meta2$Set == "D1_N_PRRSV")
#PRRSV = 10
sum(meta2$Set == "D1_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D1.pc = results(SRD129.D1.nw.De, contrast=c("Set","D1_N_PRRSV","D1_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D1.pc
sigtab.D1.pc = res.D1.pc[which(res.D1.pc$padj < .05), ]
sigtab.D1.pc = cbind(as(sigtab.D1.pc, "data.frame"), as(tax_table(SRD129.D1.nw)[rownames(sigtab.D1.pc), ], "matrix"))
format(sigtab.D1.pc$padj, scientific = TRUE)
sigtab.D1.pc$newp <- format(round(sigtab.D1.pc$padj, digits = 3), scientific = TRUE)
sigtab.D1.pc$Treatment <- ifelse(sigtab.D1.pc$log2FoldChange >=0, "PRRSV", "Control")

#Summarize sigtab.D1.pc
sum.sigtab.D1.pc <- summary(sigtab.D1.pc)
sum.sigtab.D1.pc

#ggplot
deseq.D1.pc <- ggplot(sigtab.D1.pc, aes(x=reorder(rownames(sigtab.D1.pc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D1.pc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in PRRS Virus Group Relative to Control at Nasal Site on Day 1')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "PRRSV"), values = c('#999999', '#99FF66'))
deseq.D1.pc

#Add OTU and comparisons columns
sigtab.D1.pc
sigtab.D1.pc$OTU <- rownames(sigtab.D1.pc)
sigtab.D1.pc
sigtab.D1.pc$comp <- 'D1_nasal_PRRSVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D1.pc)







################################################## Day 3 ############################################################

############## Day 3 Nasal #########################

sample_data(SRD129.genus)
SRD129.D3.nw <- subset_samples(SRD129.genus, Day == 'D3' & Tissue == 'N')
sample_sums(SRD129.D3.nw)
colnames(otu_table(SRD129.D3.nw)) #check on all the sample names
SRD129.D3.nw <- prune_taxa(taxa_sums(SRD129.D3.nw) > 1, SRD129.D3.nw)
#if taxa_sums is >1, then it will print that out in SRD129.D3.nw object and not include anything with <1.
rowSums(SRD129.D3.nw@otu_table)
SRD129.D3.nw.De <- phyloseq_to_deseq2(SRD129.D3.nw, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
SRD129.D3.nw.De <- DESeq(SRD129.D3.nw.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing


######### Day 3 Nasal BB vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D3_N_BB")
#BB = 9
sum(meta2$Set == "D3_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D3.bc = results(SRD129.D3.nw.De, contrast=c("Set","D3_N_BB","D3_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D3.bc
sigtab.D3.bc = res.D3.bc[which(res.D3.bc$padj < .05), ]
sigtab.D3.bc = cbind(as(sigtab.D3.bc, "data.frame"), as(tax_table(SRD129.D3.nw)[rownames(sigtab.D3.bc), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D3.bc$padj, scientific = TRUE)
sigtab.D3.bc$newp <- format(round(sigtab.D3.bc$padj, digits = 3), scientific = TRUE)
sigtab.D3.bc$Treatment <- ifelse(sigtab.D3.bc$log2FoldChange >=0, "BB", "Control")
#ifelse: Assigning "In-feed" = yes, "Control" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D3.bc

#Summarize sigtab.D3.bc
sum.sigtab.D3.bc <- summary(sigtab.D3.bc)
sum.sigtab.D3.bc

#ggplot
deseq.D3.bc <- ggplot(sigtab.D3.bc, aes(x=reorder(rownames(sigtab.D3.bc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D3.bc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group Relative to Control at Nasal Site on Day 3')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("BB", "Control"), values = c('#E69F00', '#999999'))
deseq.D3.bc

#Add OTU and comparisons columns
sigtab.D3.bc$OTU <- rownames(sigtab.D3.bc)
sigtab.D3.bc$comp <- 'D3_nasal_BBvsControl'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D3.bc)


######### Day 3 Nasal IAV vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D3_N_IAV")
#IAV = 10
sum(meta2$Set == "D3_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D3.nw.De$Set
res.D3.ic = results(SRD129.D3.nw.De, contrast=c("Set","D3_N_IAV","D3_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(SRD129.D3.nw.De, contrast=c("Set","D3_N_IAV","D3_N_Control")) 
sigtab.D3.ic = res.D3.ic[which(res.D3.ic$padj < .05), ]
sigtab.D3.ic = cbind(as(sigtab.D3.ic, "data.frame"), as(tax_table(SRD129.D3.nw)[rownames(sigtab.D3.ic), ], "matrix"))
format(sigtab.D3.ic$padj, scientific = TRUE)
sigtab.D3.ic$newp <- format(round(sigtab.D3.ic$padj, digits = 3), scientific = TRUE)
sigtab.D3.ic$Treatment <- ifelse(sigtab.D3.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D3.ic
sum.sigtab.D3.ic <- summary(sigtab.D3.ic)
sum.sigtab.D3.ic

#ggplot
deseq.D3.ic <- ggplot(sigtab.D3.ic, aes(x=reorder(rownames(sigtab.D3.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D3.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 3')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "IAV"), values = c('#999999', '#CC0066'))
deseq.D3.ic

#Add OTU and comparisons columns
sigtab.D3.ic
sigtab.D3.ic$OTU <- rownames(sigtab.D3.ic)
sigtab.D3.ic
sigtab.D3.ic$comp <- 'D3_nasal_IAVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D3.ic)


######### Day 3 Nasal PRRSV vs Control ###################

meta2$Set
#Number of pigs per group: 
sum(meta2$Set == "D3_N_PRRSV")
#PRRSV = 10
sum(meta2$Set == "D3_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D3.pc = results(SRD129.D3.nw.De, contrast=c("Set","D3_N_PRRSV","D3_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D3.pc
sigtab.D3.pc = res.D3.pc[which(res.D3.pc$padj < .05), ]
sigtab.D3.pc = cbind(as(sigtab.D3.pc, "data.frame"), as(tax_table(SRD129.D3.nw)[rownames(sigtab.D3.pc), ], "matrix"))
format(sigtab.D3.pc$padj, scientific = TRUE)
sigtab.D3.pc$newp <- format(round(sigtab.D3.pc$padj, digits = 3), scientific = TRUE)
sigtab.D3.pc$Treatment <- ifelse(sigtab.D3.pc$log2FoldChange >=0, "PRRSV", "Control")

#Summarize sigtab.D3.pc
sum.sigtab.D3.pc <- summary(sigtab.D3.pc)
sum.sigtab.D3.pc

#ggplot
deseq.D3.pc <- ggplot(sigtab.D3.pc, aes(x=reorder(rownames(sigtab.D3.pc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D3.pc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in PRRS Virus Group Relative to Control at Nasal Site on Day 3')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "PRRSV"), values = c('#999999', '#99FF66'))
deseq.D3.pc

#Add OTU and comparisons columns
sigtab.D3.pc
sigtab.D3.pc$OTU <- rownames(sigtab.D3.pc)
sigtab.D3.pc
sigtab.D3.pc$comp <- 'D3_nasal_PRRSVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D3.pc)




################################################## Day 7 ############################################################

############## Day 7 Nasal #########################

sample_data(SRD129.genus)
SRD129.D7.nw <- subset_samples(SRD129.genus, Day == 'D7' & Tissue == 'N')
sample_sums(SRD129.D7.nw)
colnames(otu_table(SRD129.D7.nw)) #check on all the sample names
SRD129.D7.nw <- prune_taxa(taxa_sums(SRD129.D7.nw) > 1, SRD129.D7.nw)
#if taxa_sums is >1, then it will print that out in SRD129.D7.nw object and not include anything with <1.
rowSums(SRD129.D7.nw@otu_table)
SRD129.D7.nw.De <- phyloseq_to_deseq2(SRD129.D7.nw, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
SRD129.D7.nw.De <- DESeq(SRD129.D7.nw.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing


######### Day 7 Nasal BB vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D7_N_BB")
#BB = 9
sum(meta2$Set == "D7_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D7.bc = results(SRD129.D7.nw.De, contrast=c("Set","D7_N_BB","D7_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D7.bc
sigtab.D7.bc = res.D7.bc[which(res.D7.bc$padj < .05), ]
sigtab.D7.bc = cbind(as(sigtab.D7.bc, "data.frame"), as(tax_table(SRD129.D7.nw)[rownames(sigtab.D7.bc), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D7.bc$padj, scientific = TRUE)
sigtab.D7.bc$newp <- format(round(sigtab.D7.bc$padj, digits = 3), scientific = TRUE)
sigtab.D7.bc$Treatment <- ifelse(sigtab.D7.bc$log2FoldChange >=0, "BB", "Control")
#ifelse: Assigning "In-feed" = yes, "Control" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D7.bc

#Summarize sigtab.D7.bc
sum.sigtab.D7.bc <- summary(sigtab.D7.bc)
sum.sigtab.D7.bc

#ggplot
deseq.D7.bc <- ggplot(sigtab.D7.bc, aes(x=reorder(rownames(sigtab.D7.bc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.bc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group Relative to Control at Nasal Site on Day 7')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("BB", "Control"), values = c('#E69F00', '#999999'))
deseq.D7.bc


#Add OTU and comparisons columns
sigtab.D7.bc$OTU <- rownames(sigtab.D7.bc)
sigtab.D7.bc$comp <- 'D7_nasal_BBvsControl'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D7.bc)


######### Day 7 Nasal IAV vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D7_N_IAV")
#IAV = 10
sum(meta2$Set == "D7_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D7.nw.De$Set
res.D7.ic = results(SRD129.D7.nw.De, contrast=c("Set","D7_N_IAV","D7_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(SRD129.D7.nw.De, contrast=c("Set","D7_N_IAV","D7_N_Control")) 
sigtab.D7.ic = res.D7.ic[which(res.D7.ic$padj < .05), ]
sigtab.D7.ic = cbind(as(sigtab.D7.ic, "data.frame"), as(tax_table(SRD129.D7.nw)[rownames(sigtab.D7.ic), ], "matrix"))
format(sigtab.D7.ic$padj, scientific = TRUE)
sigtab.D7.ic$newp <- format(round(sigtab.D7.ic$padj, digits = 3), scientific = TRUE)
sigtab.D7.ic$Treatment <- ifelse(sigtab.D7.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D7.ic
sum.sigtab.D7.ic <- summary(sigtab.D7.ic)
sum.sigtab.D7.ic

#ggplot
deseq.D7.ic <- ggplot(sigtab.D7.ic, aes(x=reorder(rownames(sigtab.D7.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "IAV"), values = c('#999999', '#CC0066'))
deseq.D7.ic

#Add OTU and comparisons columns
sigtab.D7.ic
sigtab.D7.ic$OTU <- rownames(sigtab.D7.ic)
sigtab.D7.ic
sigtab.D7.ic$comp <- 'D7_nasal_IAVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D7.ic)


######### Day 7 Nasal PRRSV vs Control ###################

meta2$Set
#Number of pigs per group: 
sum(meta2$Set == "D7_N_PRRSV")
#PRRSV = 10
sum(meta2$Set == "D7_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D7.pc = results(SRD129.D7.nw.De, contrast=c("Set","D7_N_PRRSV","D7_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D7.pc
sigtab.D7.pc = res.D7.pc[which(res.D7.pc$padj < .05), ]
sigtab.D7.pc = cbind(as(sigtab.D7.pc, "data.frame"), as(tax_table(SRD129.D7.nw)[rownames(sigtab.D7.pc), ], "matrix"))
format(sigtab.D7.pc$padj, scientific = TRUE)
sigtab.D7.pc$newp <- format(round(sigtab.D7.pc$padj, digits = 3), scientific = TRUE)
sigtab.D7.pc$Treatment <- ifelse(sigtab.D7.pc$log2FoldChange >=0, "PRRSV", "Control")

#Summarize sigtab.D7.pc
sum.sigtab.D7.pc <- summary(sigtab.D7.pc)
sum.sigtab.D7.pc

#ggplot
deseq.D7.pc <- ggplot(sigtab.D7.pc, aes(x=reorder(rownames(sigtab.D7.pc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.pc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in PRRS Virus Group Relative to Control at Nasal Site on Day 7')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "PRRSV"), values = c('#999999', '#99FF66'))
deseq.D7.pc

#Add OTU and comparisons columns
sigtab.D7.pc
sigtab.D7.pc$OTU <- rownames(sigtab.D7.pc)
sigtab.D7.pc
sigtab.D7.pc$comp <- 'D7_nasal_PRRSVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D7.pc)





################################################## Day 10 ############################################################

############## Day 10 Nasal #########################

sample_data(SRD129.genus)
SRD129.D10.nw <- subset_samples(SRD129.genus, Day == 'D10' & Tissue == 'N')
sample_sums(SRD129.D10.nw)
colnames(otu_table(SRD129.D10.nw)) #check on all the sample names
SRD129.D10.nw <- prune_taxa(taxa_sums(SRD129.D10.nw) > 1, SRD129.D10.nw)
#if taxa_sums is >1, then it will print that out in SRD129.D10.nw object and not include anything with <1.
rowSums(SRD129.D10.nw@otu_table)
SRD129.D10.nw.De <- phyloseq_to_deseq2(SRD129.D10.nw, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
SRD129.D10.nw.De <- DESeq(SRD129.D10.nw.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing


######### Day 10 Nasal BB vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D10_N_BB")
#BB = 10
sum(meta2$Set == "D10_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D10.bc = results(SRD129.D10.nw.De, contrast=c("Set","D10_N_BB","D10_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D10.bc
sigtab.D10.bc = res.D10.bc[which(res.D10.bc$padj < .05), ]
sigtab.D10.bc = cbind(as(sigtab.D10.bc, "data.frame"), as(tax_table(SRD129.D10.nw)[rownames(sigtab.D10.bc), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D10.bc$padj, scientific = TRUE)
sigtab.D10.bc$newp <- format(round(sigtab.D10.bc$padj, digits = 3), scientific = TRUE)
sigtab.D10.bc$Treatment <- ifelse(sigtab.D10.bc$log2FoldChange >=0, "BB", "Control")
#ifelse: Assigning "In-feed" = yes, "Control" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D10.bc

#Summarize sigtab.D10.bc
sum.sigtab.D10.bc <- summary(sigtab.D10.bc)
sum.sigtab.D10.bc

#ggplot
deseq.D10.bc <- ggplot(sigtab.D10.bc, aes(x=reorder(rownames(sigtab.D10.bc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D10.bc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group Relative to Control at Nasal Site on Day 10')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("BB", "Control"), values = c('#E69F00', '#999999'))
deseq.D10.bc

#Add OTU and comparisons columns
sigtab.D10.bc$OTU <- rownames(sigtab.D10.bc)
sigtab.D10.bc$comp <- 'D10_nasal_BBvsControl'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D10.bc)


######### Day 10 Nasal IAV vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D10_N_IAV")
#IAV = 10
sum(meta2$Set == "D10_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D10.nw.De$Set
res.D10.ic = results(SRD129.D10.nw.De, contrast=c("Set","D10_N_IAV","D10_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(SRD129.D10.nw.De, contrast=c("Set","D10_N_IAV","D10_N_Control")) 
sigtab.D10.ic = res.D10.ic[which(res.D10.ic$padj < .05), ]
sigtab.D10.ic = cbind(as(sigtab.D10.ic, "data.frame"), as(tax_table(SRD129.D10.nw)[rownames(sigtab.D10.ic), ], "matrix"))
format(sigtab.D10.ic$padj, scientific = TRUE)
sigtab.D10.ic$newp <- format(round(sigtab.D10.ic$padj, digits = 3), scientific = TRUE)
sigtab.D10.ic$Treatment <- ifelse(sigtab.D10.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D10.ic
sum.sigtab.D10.ic <- summary(sigtab.D10.ic)
sum.sigtab.D10.ic

#ggplot
deseq.D10.ic <- ggplot(sigtab.D10.ic, aes(x=reorder(rownames(sigtab.D10.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D10.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 10')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "IAV"), values = c('#999999', '#CC0066'))
deseq.D10.ic

#Add OTU and comparisons columns
sigtab.D10.ic
sigtab.D10.ic$OTU <- rownames(sigtab.D10.ic)
sigtab.D10.ic
sigtab.D10.ic$comp <- 'D10_nasal_IAVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D10.ic)


######### Day 10 Nasal PRRSV vs Control ###################

meta2$Set
#Number of pigs per group: 
sum(meta2$Set == "D10_N_PRRSV")
#PRRSV = 10
sum(meta2$Set == "D10_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D10.pc = results(SRD129.D10.nw.De, contrast=c("Set","D10_N_PRRSV","D10_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D10.pc
sigtab.D10.pc = res.D10.pc[which(res.D10.pc$padj < .05), ]
sigtab.D10.pc = cbind(as(sigtab.D10.pc, "data.frame"), as(tax_table(SRD129.D10.nw)[rownames(sigtab.D10.pc), ], "matrix"))
format(sigtab.D10.pc$padj, scientific = TRUE)
sigtab.D10.pc$newp <- format(round(sigtab.D10.pc$padj, digits = 3), scientific = TRUE)
sigtab.D10.pc$Treatment <- ifelse(sigtab.D10.pc$log2FoldChange >=0, "PRRSV", "Control")

#Summarize sigtab.D10.pc
sum.sigtab.D10.pc <- summary(sigtab.D10.pc)
sum.sigtab.D10.pc

#ggplot
deseq.D10.pc <- ggplot(sigtab.D10.pc, aes(x=reorder(rownames(sigtab.D10.pc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D10.pc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in PRRS Virus Group Relative to Control at Nasal Site on Day 10')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "PRRSV"), values = c('#999999', '#99FF66'))
deseq.D10.pc

#Add OTU and comparisons columns
sigtab.D10.pc
sigtab.D10.pc$OTU <- rownames(sigtab.D10.pc)
sigtab.D10.pc
sigtab.D10.pc$comp <- 'D10_nasal_PRRSVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D10.pc)



################################################## Day 14 ############################################################

############## Day 14 Nasal #########################

sample_data(SRD129.genus)
SRD129.D14.nw <- subset_samples(SRD129.genus, Day == 'D14' & Tissue == 'N')
sample_sums(SRD129.D14.nw)
colnames(otu_table(SRD129.D14.nw)) #check on all the sample names
SRD129.D14.nw <- prune_taxa(taxa_sums(SRD129.D14.nw) > 1, SRD129.D14.nw)
#if taxa_sums is >1, then it will print that out in SRD129.D14.nw object and not include anything with <1.
rowSums(SRD129.D14.nw@otu_table)
SRD129.D14.nw.De <- phyloseq_to_deseq2(SRD129.D14.nw, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
SRD129.D14.nw.De <- DESeq(SRD129.D14.nw.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing


######### Day 14 Nasal BB vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D14_N_BB")
#BB = 10
sum(meta2$Set == "D14_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D14.bc = results(SRD129.D14.nw.De, contrast=c("Set","D14_N_BB","D14_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D14.bc
sigtab.D14.bc = res.D14.bc[which(res.D14.bc$padj < .05), ]
sigtab.D14.bc = cbind(as(sigtab.D14.bc, "data.frame"), as(tax_table(SRD129.D14.nw)[rownames(sigtab.D14.bc), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D14.bc$padj, scientific = TRUE)
sigtab.D14.bc$newp <- format(round(sigtab.D14.bc$padj, digits = 3), scientific = TRUE)
sigtab.D14.bc$Treatment <- ifelse(sigtab.D14.bc$log2FoldChange >=0, "BB", "Control")
#ifelse: Assigning "In-feed" = yes, "Control" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D14.bc

#Summarize sigtab.D14.bc
sum.sigtab.D14.bc <- summary(sigtab.D14.bc)
sum.sigtab.D14.bc

#ggplot
deseq.D14.bc <- ggplot(sigtab.D14.bc, aes(x=reorder(rownames(sigtab.D14.bc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.bc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group Relative to Control at Nasal Site on Day 14')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("BB", "Control"), values = c('#E69F00', '#999999'))
deseq.D14.bc

#Add OTU and comparisons columns
sigtab.D14.bc$OTU <- rownames(sigtab.D14.bc)
sigtab.D14.bc$comp <- 'D14_nasal_BBvsControl'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D14.bc)


######### Day 14 Nasal IAV vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D14_N_IAV")
#IAV = 10
sum(meta2$Set == "D14_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D14.nw.De$Set
res.D14.ic = results(SRD129.D14.nw.De, contrast=c("Set","D14_N_IAV","D14_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(SRD129.D14.nw.De, contrast=c("Set","D14_N_IAV","D14_N_Control")) 
sigtab.D14.ic = res.D14.ic[which(res.D14.ic$padj < .05), ]
sigtab.D14.ic = cbind(as(sigtab.D14.ic, "data.frame"), as(tax_table(SRD129.D14.nw)[rownames(sigtab.D14.ic), ], "matrix"))
format(sigtab.D14.ic$padj, scientific = TRUE)
sigtab.D14.ic$newp <- format(round(sigtab.D14.ic$padj, digits = 3), scientific = TRUE)
sigtab.D14.ic$Treatment <- ifelse(sigtab.D14.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D14.ic
sum.sigtab.D14.ic <- summary(sigtab.D14.ic)
sum.sigtab.D14.ic

#ggplot
deseq.D14.ic <- ggplot(sigtab.D14.ic, aes(x=reorder(rownames(sigtab.D14.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 14')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "IAV"), values = c('#999999', '#CC0066'))
deseq.D14.ic

#Add OTU and comparisons columns
sigtab.D14.ic
sigtab.D14.ic$OTU <- rownames(sigtab.D14.ic)
sigtab.D14.ic
sigtab.D14.ic$comp <- 'D14_nasal_IAVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D14.ic)


######### Day 14 Nasal PRRSV vs Control ###################

meta2$Set
#Number of pigs per group: 
sum(meta2$Set == "D14_N_PRRSV")
#PRRSV = 10
sum(meta2$Set == "D14_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D14.pc = results(SRD129.D14.nw.De, contrast=c("Set","D14_N_PRRSV","D14_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D14.pc
sigtab.D14.pc = res.D14.pc[which(res.D14.pc$padj < .05), ]
sigtab.D14.pc = cbind(as(sigtab.D14.pc, "data.frame"), as(tax_table(SRD129.D14.nw)[rownames(sigtab.D14.pc), ], "matrix"))
format(sigtab.D14.pc$padj, scientific = TRUE)
sigtab.D14.pc$newp <- format(round(sigtab.D14.pc$padj, digits = 3), scientific = TRUE)
sigtab.D14.pc$Treatment <- ifelse(sigtab.D14.pc$log2FoldChange >=0, "PRRSV", "Control")

#Summarize sigtab.D14.pc
sum.sigtab.D14.pc <- summary(sigtab.D14.pc)
sum.sigtab.D14.pc

#ggplot
deseq.D14.pc <- ggplot(sigtab.D14.pc, aes(x=reorder(rownames(sigtab.D14.pc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.pc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in PRRS Virus Group Relative to Control at Nasal Site on Day 14')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "PRRSV"), values = c('#999999', '#99FF66'))
deseq.D14.pc

#Add OTU and comparisons columns
sigtab.D14.pc
sigtab.D14.pc$OTU <- rownames(sigtab.D14.pc)
sigtab.D14.pc
sigtab.D14.pc$comp <- 'D14_nasal_PRRSVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D14.pc)



################################################## Day 21 ############################################################

############## Day 21 Nasal #########################

sample_data(SRD129.genus)
SRD129.D21.nw <- subset_samples(SRD129.genus, Day == 'D21' & Tissue == 'N')
sample_sums(SRD129.D21.nw)
colnames(otu_table(SRD129.D21.nw)) #check on all the sample names
SRD129.D21.nw <- prune_taxa(taxa_sums(SRD129.D21.nw) > 1, SRD129.D21.nw)
#if taxa_sums is >1, then it will print that out in SRD129.D21.nw object and not include anything with <1.
rowSums(SRD129.D21.nw@otu_table)
SRD129.D21.nw.De <- phyloseq_to_deseq2(SRD129.D21.nw, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
SRD129.D21.nw.De <- DESeq(SRD129.D21.nw.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing


######### Day 21 Nasal BB vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D21_N_BB")
#BB = 10
sum(meta2$Set == "D21_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D21.bc = results(SRD129.D21.nw.De, contrast=c("Set","D21_N_BB","D21_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D21.bc
sigtab.D21.bc = res.D21.bc[which(res.D21.bc$padj < .05), ]
sigtab.D21.bc = cbind(as(sigtab.D21.bc, "data.frame"), as(tax_table(SRD129.D21.nw)[rownames(sigtab.D21.bc), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D21.bc$padj, scientific = TRUE)
sigtab.D21.bc$newp <- format(round(sigtab.D21.bc$padj, digits = 3), scientific = TRUE)
sigtab.D21.bc$Treatment <- ifelse(sigtab.D21.bc$log2FoldChange >=0, "BB", "Control")
#ifelse: Assigning "In-feed" = yes, "Control" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D21.bc

#Summarize sigtab.D21.bc
sum.sigtab.D21.bc <- summary(sigtab.D21.bc)
sum.sigtab.D21.bc

#ggplot
deseq.D21.bc <- ggplot(sigtab.D21.bc, aes(x=reorder(rownames(sigtab.D21.bc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D21.bc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group Relative to Control at Nasal Site on Day 21')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("BB", "Control"), values = c('#E69F00', '#999999'))
deseq.D21.bc

#Add OTU and comparisons columns
sigtab.D21.bc$OTU <- rownames(sigtab.D21.bc)
sigtab.D21.bc$comp <- 'D21_nasal_BBvsControl'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D21.bc)


######### Day 21 Nasal IAV vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D21_N_IAV")
#IAV = 9
sum(meta2$Set == "D21_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D21.nw.De$Set
res.D21.ic = results(SRD129.D21.nw.De, contrast=c("Set","D21_N_IAV","D21_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(SRD129.D21.nw.De, contrast=c("Set","D21_N_IAV","D21_N_Control")) 
sigtab.D21.ic = res.D21.ic[which(res.D21.ic$padj < .05), ]
sigtab.D21.ic = cbind(as(sigtab.D21.ic, "data.frame"), as(tax_table(SRD129.D21.nw)[rownames(sigtab.D21.ic), ], "matrix"))
format(sigtab.D21.ic$padj, scientific = TRUE)
sigtab.D21.ic$newp <- format(round(sigtab.D21.ic$padj, digits = 3), scientific = TRUE)
sigtab.D21.ic$Treatment <- ifelse(sigtab.D21.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D21.ic
sum.sigtab.D21.ic <- summary(sigtab.D21.ic)
sum.sigtab.D21.ic

#ggplot
deseq.D21.ic <- ggplot(sigtab.D21.ic, aes(x=reorder(rownames(sigtab.D21.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D21.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 21')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "IAV"), values = c('#999999', '#CC0066'))
deseq.D21.ic

#Add OTU and comparisons columns
sigtab.D21.ic
sigtab.D21.ic$OTU <- rownames(sigtab.D21.ic)
sigtab.D21.ic
sigtab.D21.ic$comp <- 'D21_nasal_IAVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D21.ic)


######### Day 21 Nasal PRRSV vs Control ###################

meta2$Set
#Number of pigs per group: 
sum(meta2$Set == "D21_N_PRRSV")
#PRRSV = 10
sum(meta2$Set == "D21_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D21.pc = results(SRD129.D21.nw.De, contrast=c("Set","D21_N_PRRSV","D21_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D21.pc
sigtab.D21.pc = res.D21.pc[which(res.D21.pc$padj < .05), ]
sigtab.D21.pc = cbind(as(sigtab.D21.pc, "data.frame"), as(tax_table(SRD129.D21.nw)[rownames(sigtab.D21.pc), ], "matrix"))
format(sigtab.D21.pc$padj, scientific = TRUE)
sigtab.D21.pc$newp <- format(round(sigtab.D21.pc$padj, digits = 3), scientific = TRUE)
sigtab.D21.pc$Treatment <- ifelse(sigtab.D21.pc$log2FoldChange >=0, "PRRSV", "Control")

#Summarize sigtab.D21.pc
sum.sigtab.D21.pc <- summary(sigtab.D21.pc)
sum.sigtab.D21.pc

#ggplot
deseq.D21.pc <- ggplot(sigtab.D21.pc, aes(x=reorder(rownames(sigtab.D21.pc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D21.pc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in PRRS Virus Group Relative to Control at Nasal Site on Day 21')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "PRRSV"), values = c('#999999', '#99FF66'))
deseq.D21.pc

#Add OTU and comparisons columns
sigtab.D21.pc
sigtab.D21.pc$OTU <- rownames(sigtab.D21.pc)
sigtab.D21.pc
sigtab.D21.pc$comp <- 'D21_nasal_PRRSVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D21.pc)







################################################## Day 28 ############################################################

############## Day 28 Nasal #########################

sample_data(SRD129.genus)
SRD129.D28.nw <- subset_samples(SRD129.genus, Day == 'D28' & Tissue == 'N')
sample_sums(SRD129.D28.nw)
colnames(otu_table(SRD129.D28.nw)) #check on all the sample names
SRD129.D28.nw <- prune_taxa(taxa_sums(SRD129.D28.nw) > 1, SRD129.D28.nw)
#if taxa_sums is >1, then it will print that out in SRD129.D28.nw object and not include anything with <1.
rowSums(SRD129.D28.nw@otu_table)
SRD129.D28.nw.De <- phyloseq_to_deseq2(SRD129.D28.nw, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
SRD129.D28.nw.De <- DESeq(SRD129.D28.nw.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing


######### Day 28 Nasal BB vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D28_N_BB")
#BB = 10
sum(meta2$Set == "D28_N_Control")
#Control = 4

#Extract results from a DESeq analysis, organize table
res.D28.bc = results(SRD129.D28.nw.De, contrast=c("Set","D28_N_BB","D28_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D28.bc
sigtab.D28.bc = res.D28.bc[which(res.D28.bc$padj < .05), ]
sigtab.D28.bc = cbind(as(sigtab.D28.bc, "data.frame"), as(tax_table(SRD129.D28.nw)[rownames(sigtab.D28.bc), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D28.bc$padj, scientific = TRUE)
sigtab.D28.bc$newp <- format(round(sigtab.D28.bc$padj, digits = 3), scientific = TRUE)
sigtab.D28.bc$Treatment <- ifelse(sigtab.D28.bc$log2FoldChange >=0, "BB", "Control")
#ifelse: Assigning "In-feed" = yes, "Control" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D28.bc

#Summarize sigtab.D28.bc
sum.sigtab.D28.bc <- summary(sigtab.D28.bc)
sum.sigtab.D28.bc

#ggplot
deseq.D28.bc <- ggplot(sigtab.D28.bc, aes(x=reorder(rownames(sigtab.D28.bc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D28.bc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group Relative to Control at Nasal Site on Day 28')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("BB", "Control"), values = c('#E69F00', '#999999'))
deseq.D28.bc

#Add OTU and comparisons columns
sigtab.D28.bc$OTU <- rownames(sigtab.D28.bc)
sigtab.D28.bc$comp <- 'D28_nasal_BBvsControl'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D28.bc)


######### Day 28 Nasal IAV vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D28_N_IAV")
#IAV = 5
sum(meta2$Set == "D28_N_Control")
#Control = 4

#Extract results from a DESeq analysis, organize table
SRD129.D28.nw.De$Set
res.D28.ic = results(SRD129.D28.nw.De, contrast=c("Set","D28_N_IAV","D28_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(SRD129.D28.nw.De, contrast=c("Set","D28_N_IAV","D28_N_Control")) 
sigtab.D28.ic = res.D28.ic[which(res.D28.ic$padj < .05), ]
sigtab.D28.ic = cbind(as(sigtab.D28.ic, "data.frame"), as(tax_table(SRD129.D28.nw)[rownames(sigtab.D28.ic), ], "matrix"))
format(sigtab.D28.ic$padj, scientific = TRUE)
sigtab.D28.ic$newp <- format(round(sigtab.D28.ic$padj, digits = 3), scientific = TRUE)
sigtab.D28.ic$Treatment <- ifelse(sigtab.D28.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D28.ic
sum.sigtab.D28.ic <- summary(sigtab.D28.ic)
sum.sigtab.D28.ic

#ggplot
deseq.D28.ic <- ggplot(sigtab.D28.ic, aes(x=reorder(rownames(sigtab.D28.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D28.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 28')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "IAV"), values = c('#999999', '#CC0066'))
deseq.D28.ic

#Add OTU and comparisons columns
sigtab.D28.ic
sigtab.D28.ic$OTU <- rownames(sigtab.D28.ic)
sigtab.D28.ic
sigtab.D28.ic$comp <- 'D28_nasal_IAVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D28.ic)


######### Day 28 Nasal PRRSV vs Control ###################

meta2$Set
#Number of pigs per group: 
sum(meta2$Set == "D28_N_PRRSV")
#PRRSV = 8
sum(meta2$Set == "D28_N_Control")
#Control = 4

#Extract results from a DESeq analysis, organize table
res.D28.pc = results(SRD129.D28.nw.De, contrast=c("Set","D28_N_PRRSV","D28_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D28.pc
sigtab.D28.pc = res.D28.pc[which(res.D28.pc$padj < .05), ]
sigtab.D28.pc = cbind(as(sigtab.D28.pc, "data.frame"), as(tax_table(SRD129.D28.nw)[rownames(sigtab.D28.pc), ], "matrix"))
format(sigtab.D28.pc$padj, scientific = TRUE)
sigtab.D28.pc$newp <- format(round(sigtab.D28.pc$padj, digits = 3), scientific = TRUE)
sigtab.D28.pc$Treatment <- ifelse(sigtab.D28.pc$log2FoldChange >=0, "PRRSV", "Control")

#Summarize sigtab.D28.pc
sum.sigtab.D28.pc <- summary(sigtab.D28.pc)
sum.sigtab.D28.pc

#ggplot
deseq.D28.pc <- ggplot(sigtab.D28.pc, aes(x=reorder(rownames(sigtab.D28.pc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D28.pc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in PRRS Virus Group Relative to Control at Nasal Site on Day 28')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "PRRSV"), values = c('#999999', '#99FF66'))
deseq.D28.pc

#Add OTU and comparisons columns
sigtab.D28.pc
sigtab.D28.pc$OTU <- rownames(sigtab.D28.pc)
sigtab.D28.pc
sigtab.D28.pc$comp <- 'D28_nasal_PRRSVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D28.pc)




################################################## Day 36 ############################################################

############## Day 36 Nasal #########################

sample_data(SRD129.genus)
SRD129.D36.nw <- subset_samples(SRD129.genus, Day == 'D36' & Tissue == 'N')
sample_sums(SRD129.D36.nw)
colnames(otu_table(SRD129.D36.nw)) #check on all the sample names
SRD129.D36.nw <- prune_taxa(taxa_sums(SRD129.D36.nw) > 1, SRD129.D36.nw)
#if taxa_sums is >1, then it will print that out in SRD129.D36.nw object and not include anything with <1.
rowSums(SRD129.D36.nw@otu_table)
SRD129.D36.nw.De <- phyloseq_to_deseq2(SRD129.D36.nw, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
SRD129.D36.nw.De <- DESeq(SRD129.D36.nw.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing


######### Day 36 Nasal BB vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D36_N_BB")
#BB = 10
sum(meta2$Set == "D36_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D36.bc = results(SRD129.D36.nw.De, contrast=c("Set","D36_N_BB","D36_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D36.bc
sigtab.D36.bc = res.D36.bc[which(res.D36.bc$padj < .05), ]
sigtab.D36.bc = cbind(as(sigtab.D36.bc, "data.frame"), as(tax_table(SRD129.D36.nw)[rownames(sigtab.D36.bc), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D36.bc$padj, scientific = TRUE)
sigtab.D36.bc$newp <- format(round(sigtab.D36.bc$padj, digits = 3), scientific = TRUE)
sigtab.D36.bc$Treatment <- ifelse(sigtab.D36.bc$log2FoldChange >=0, "BB", "Control")
#ifelse: Assigning "In-feed" = yes, "Control" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D36.bc

#Summarize sigtab.D36.bc
sum.sigtab.D36.bc <- summary(sigtab.D36.bc)
sum.sigtab.D36.bc

#ggplot
deseq.D36.bc <- ggplot(sigtab.D36.bc, aes(x=reorder(rownames(sigtab.D36.bc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D36.bc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group Relative to Control at Nasal Site on Day 36')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("BB", "Control"), values = c('#E69F00', '#999999'))
deseq.D36.bc

#Add OTU and comparisons columns
sigtab.D36.bc$OTU <- rownames(sigtab.D36.bc)
sigtab.D36.bc$comp <- 'D36_nasal_BBvsControl'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D36.bc)


######### Day 36 Nasal IAV vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D36_N_IAV")
#IAV = 9
sum(meta2$Set == "D36_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D36.nw.De$Set
res.D36.ic = results(SRD129.D36.nw.De, contrast=c("Set","D36_N_IAV","D36_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(SRD129.D36.nw.De, contrast=c("Set","D36_N_IAV","D36_N_Control")) 
sigtab.D36.ic = res.D36.ic[which(res.D36.ic$padj < .05), ]
sigtab.D36.ic = cbind(as(sigtab.D36.ic, "data.frame"), as(tax_table(SRD129.D36.nw)[rownames(sigtab.D36.ic), ], "matrix"))
format(sigtab.D36.ic$padj, scientific = TRUE)
sigtab.D36.ic$newp <- format(round(sigtab.D36.ic$padj, digits = 3), scientific = TRUE)
sigtab.D36.ic$Treatment <- ifelse(sigtab.D36.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D36.ic
sum.sigtab.D36.ic <- summary(sigtab.D36.ic)
sum.sigtab.D36.ic

#ggplot
deseq.D36.ic <- ggplot(sigtab.D36.ic, aes(x=reorder(rownames(sigtab.D36.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D36.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 36')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "IAV"), values = c('#999999', '#CC0066'))
deseq.D36.ic

#Add OTU and comparisons columns
sigtab.D36.ic
sigtab.D36.ic$OTU <- rownames(sigtab.D36.ic)
sigtab.D36.ic
sigtab.D36.ic$comp <- 'D36_nasal_IAVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D36.ic)


######### Day 36 Nasal PRRSV vs Control ###################

meta2$Set
#Number of pigs per group: 
sum(meta2$Set == "D36_N_PRRSV")
#PRRSV = 10
sum(meta2$Set == "D36_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D36.pc = results(SRD129.D36.nw.De, contrast=c("Set","D36_N_PRRSV","D36_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D36.pc
sigtab.D36.pc = res.D36.pc[which(res.D36.pc$padj < .05), ]
sigtab.D36.pc = cbind(as(sigtab.D36.pc, "data.frame"), as(tax_table(SRD129.D36.nw)[rownames(sigtab.D36.pc), ], "matrix"))
format(sigtab.D36.pc$padj, scientific = TRUE)
sigtab.D36.pc$newp <- format(round(sigtab.D36.pc$padj, digits = 3), scientific = TRUE)
sigtab.D36.pc$Treatment <- ifelse(sigtab.D36.pc$log2FoldChange >=0, "PRRSV", "Control")

#Summarize sigtab.D36.pc
sum.sigtab.D36.pc <- summary(sigtab.D36.pc)
sum.sigtab.D36.pc

#ggplot
deseq.D36.pc <- ggplot(sigtab.D36.pc, aes(x=reorder(rownames(sigtab.D36.pc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D36.pc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in PRRS Virus Group Relative to Control at Nasal Site on Day 36')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "PRRSV"), values = c('#999999', '#99FF66'))
deseq.D36.pc

#Add OTU and comparisons columns
sigtab.D36.pc
sigtab.D36.pc$OTU <- rownames(sigtab.D36.pc)
sigtab.D36.pc
sigtab.D36.pc$comp <- 'D36_nasal_PRRSVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D36.pc)


################################################## Day 42 ############################################################

############## Day 42 Nasal #########################

sample_data(SRD129.genus)
SRD129.D42.nw <- subset_samples(SRD129.genus, Day == 'D42' & Tissue == 'N')
sample_sums(SRD129.D42.nw)
colnames(otu_table(SRD129.D42.nw)) #check on all the sample names
SRD129.D42.nw <- prune_taxa(taxa_sums(SRD129.D42.nw) > 1, SRD129.D42.nw)
#if taxa_sums is >1, then it will print that out in SRD129.D42.nw object and not include anything with <1.
rowSums(SRD129.D42.nw@otu_table)
SRD129.D42.nw.De <- phyloseq_to_deseq2(SRD129.D42.nw, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
SRD129.D42.nw.De <- DESeq(SRD129.D42.nw.De, test = "Wald", fitType = "parametric")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 7 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing


######### Day 42 Nasal BB vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D42_N_BB")
#BB = 10
sum(meta2$Set == "D42_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D42.bc = results(SRD129.D42.nw.De, contrast=c("Set","D42_N_BB","D42_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D42.bc
sigtab.D42.bc = res.D42.bc[which(res.D42.bc$padj < .05), ]
sigtab.D42.bc = cbind(as(sigtab.D42.bc, "data.frame"), as(tax_table(SRD129.D42.nw)[rownames(sigtab.D42.bc), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D42.bc$padj, scientific = TRUE)
sigtab.D42.bc$newp <- format(round(sigtab.D42.bc$padj, digits = 3), scientific = TRUE)
sigtab.D42.bc$Treatment <- ifelse(sigtab.D42.bc$log2FoldChange >=0, "BB", "Control")
#ifelse: Assigning "In-feed" = yes, "Control" = no; it's important to make sure you have the
#correct group names in the "yes" and "no" position in ifelse function
sigtab.D42.bc

#Summarize sigtab.D42.bc
sum.sigtab.D42.bc <- summary(sigtab.D42.bc)
sum.sigtab.D42.bc

#ggplot
deseq.D42.bc <- ggplot(sigtab.D42.bc, aes(x=reorder(rownames(sigtab.D42.bc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D42.bc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Bordetella bronchiseptica Group Relative to Control at Nasal Site on Day 42')+ coord_flip() + 
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("BB", "Control"), values = c('#E69F00', '#999999'))
deseq.D42.bc

#Add OTU and comparisons columns
sigtab.D42.bc$OTU <- rownames(sigtab.D42.bc)
sigtab.D42.bc$comp <- 'D42_nasal_BBvsControl'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D42.bc)


######### Day 42 Nasal IAV vs Control ###################

meta2$Set
#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D42_N_IAV")
#IAV = 9
sum(meta2$Set == "D42_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
SRD129.D42.nw.De$Set
res.D42.ic = results(SRD129.D42.nw.De, contrast=c("Set","D42_N_IAV","D42_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(SRD129.D42.nw.De, contrast=c("Set","D42_N_IAV","D42_N_Control")) 
sigtab.D42.ic = res.D42.ic[which(res.D42.ic$padj < .05), ]
sigtab.D42.ic = cbind(as(sigtab.D42.ic, "data.frame"), as(tax_table(SRD129.D42.nw)[rownames(sigtab.D42.ic), ], "matrix"))
format(sigtab.D42.ic$padj, scientific = TRUE)
sigtab.D42.ic$newp <- format(round(sigtab.D42.ic$padj, digits = 3), scientific = TRUE)
sigtab.D42.ic$Treatment <- ifelse(sigtab.D42.ic$log2FoldChange >=0, "IAV", "Control")

#Summarize sigtab.D42.ic
sum.sigtab.D42.ic <- summary(sigtab.D42.ic)
sum.sigtab.D42.ic

#ggplot
deseq.D42.ic <- ggplot(sigtab.D42.ic, aes(x=reorder(rownames(sigtab.D42.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D42.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in Influenza A Virus Group Relative to Control at Nasal Site on Day 42')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "IAV"), values = c('#999999', '#CC0066'))
deseq.D42.ic

#Add OTU and comparisons columns
sigtab.D42.ic
sigtab.D42.ic$OTU <- rownames(sigtab.D42.ic)
sigtab.D42.ic
sigtab.D42.ic$comp <- 'D42_nasal_IAVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D42.ic)


######### Day 42 Nasal PRRSV vs Control ###################

meta2$Set
#Number of pigs per group: 
sum(meta2$Set == "D42_N_PRRSV")
#PRRSV = 10
sum(meta2$Set == "D42_N_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D42.pc = results(SRD129.D42.nw.De, contrast=c("Set","D42_N_PRRSV","D42_N_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D42.pc
sigtab.D42.pc = res.D42.pc[which(res.D42.pc$padj < .05), ]
sigtab.D42.pc = cbind(as(sigtab.D42.pc, "data.frame"), as(tax_table(SRD129.D42.nw)[rownames(sigtab.D42.pc), ], "matrix"))
format(sigtab.D42.pc$padj, scientific = TRUE)
sigtab.D42.pc$newp <- format(round(sigtab.D42.pc$padj, digits = 3), scientific = TRUE)
sigtab.D42.pc$Treatment <- ifelse(sigtab.D42.pc$log2FoldChange >=0, "PRRSV", "Control")

#Summarize sigtab.D42.pc
sum.sigtab.D42.pc <- summary(sigtab.D42.pc)
sum.sigtab.D42.pc

#ggplot
deseq.D42.pc <- ggplot(sigtab.D42.pc, aes(x=reorder(rownames(sigtab.D42.pc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D42.pc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
  theme(axis.text.x=element_text(color = 'black', size = 13),
        axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in PRRS Virus Group Relative to Control at Nasal Site on Day 42')+ coord_flip() +
  theme(plot.title = element_text(size = 18, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
  scale_fill_manual(labels = c("Control", "PRRSV"), values = c('#999999', '#99FF66'))
deseq.D42.pc

#Add OTU and comparisons columns
sigtab.D42.pc
sigtab.D42.pc$OTU <- rownames(sigtab.D42.pc)
sigtab.D42.pc
sigtab.D42.pc$comp <- 'D42_nasal_PRRSVvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D42.pc)



#######################################################################################################

#write csv
write.csv(final.nonsigtab, file= "SRD129_FinalDiffAbundNasalGenus_SignificantDays.csv")

#######################################################################################################



######### Plots of Diff Abund Nasal Families and Genera Combined for each pairwise comparison ########
library("ggsci")

#BB and Control
final_bc <- sigtab.D0.bc
final_bc <- rbind(final_bc, sigtab.DNEG12.bc, sigtab.DNEG6.bc, sigtab.D1.bc, sigtab.D3.bc, sigtab.D7.bc, sigtab.D10.bc, 
                  sigtab.D14.bc, sigtab.D21.bc, sigtab.D28.bc, sigtab.D36.bc, sigtab.D42.bc)
final_bc$Family_Genus <- paste(final_bc$Family, final_bc$Genus) #create new column with Family_Genus
final_bc$comp
final_bc$comp <- factor(final_bc$comp, levels=c("DNEG12_nasal_BBvsControl", "DNEG6_nasal_BBvsControl", "D0_nasal_BBvsControl",
                                           "D1_nasal_BBvsControl", "D3_nasal_BBvsControl", "D7_nasal_BBvsControl",
                                           "D10_nasal_BBvsControl", "D14_nasal_BBvsControl", "D21_nasal_BBvsControl",
                                           "D28_nasal_BBvsControl", "D36_nasal_BBvsControl", "D42_nasal_BBvsControl"))
levels(final_bc$comp)
bc_plot <- ggplot(final_bc, aes(x=Family_Genus, log2FoldChange, fill = comp)) +
  geom_bar(stat='identity') +
  labs(x="Family Genus", y = "Total log2 Fold Change") +
  theme(axis.text.x=element_text(color = 'black', size = 18),
        axis.text.y=element_text(color = 'black', size=15), 
        axis.title.x=element_text(size = 20),
        axis.title.y=element_text(size = 20))+ 
  coord_flip() +
  scale_fill_igv(name = "comp") +
  ggtitle('Differentially Abundant Nasal Families and Genera between Bordetella bronchiseptica and Control Groups') + 
  theme(plot.title = element_text(size = 20), legend.text = element_text(size=20), legend.title = element_text(size=20))
bc_plot <- bc_plot + guides(fill=guide_legend(title="Day and treatment group"))
bc_plot

#BB and control Log2fold and Basemean plots
write.csv(final_bc, file= "SRD129_FinalDiffAbundNasalGenus_BBControl.csv")
bc <- read.csv('SRD129_FinalDiffAbundNasalGenus_BBControl.csv', header = TRUE, sep = ",")
head(bc[,1:10])
colnames(bc)
bc$DayComp <- sub('_[A-Za-z]+', '\\2', bc$comp)
unique(bc$DayComp)
bc$Day <- sub('_[A-Za-z]+', '\\2', bc$DayComp)
unique(bc$Day)
bc <- bc[!grepl("DNEG12", bc$Day),]
bc <- bc[!grepl("DNEG6", bc$Day),]
unique(bc$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
bc$Day = factor(bc$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
levels(bc$Day)


(logfoldplot <- ggplot(data=bc, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 2)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
logfoldplot

(basemeanplot <- ggplot(data=bc, aes(x=Day, y=baseMean, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, scales = "free") + ylab('baseMean') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 2)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
basemeanplot


#BB and control select genera Log2fold plots
bbc <- read.csv('SRD129_FinalDiffAbundNasalGenus_BBControl_CondensedListFinal.csv', header = TRUE, sep = ",")
head(bbc[,1:10])
colnames(bbc)
bbc$DayComp <- sub('_[A-Za-z]+', '\\2', bbc$comp)
unique(bbc$DayComp) #""D0_BBvsControl"  "D1_BBvsControl"  "D3_BBvsControl"  "D7_BBvsControl"  "D10_BBvsControl" "D14_BBvsControl" "D21_BBvsControl" "D28_BBvsControl" "D36_BBvsControl" "D42_BBvsControl"
bbc$Day <- sub('_[A-Za-z]+', '\\2', bbc$DayComp)
unique(bbc$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
bbc$Day = factor(bbc$Day, levels=c("D1","D3","D7","D10","D14","D21","D28","D36","D42"))
levels(bbc$Day) #"D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
(bbc_logfoldplot <- ggplot(data=bbc, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 7, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
bbc_logfoldplot



#BB and control even smaller select genera Log2fold plots
bbd <- read.csv('SRD129_FinalDiffAbundNasalGenus_BBControl_CondensedSmallerList.csv', header = TRUE, sep = ",")
head(bbd[,1:10])
colnames(bbd)
bbd$DayComp <- sub('_[A-Za-z]+', '\\2', bbd$comp)
unique(bbd$DayComp) #""D0_BBvsControl"  "D1_BBvsControl"  "D3_BBvsControl"  "D7_BBvsControl"  "D10_BBvsControl" "D14_BBvsControl" "D21_BBvsControl" "D28_BBvsControl" "D36_BBvsControl" "D42_BBvsControl"
bbd$Day <- sub('_[A-Za-z]+', '\\2', bbd$DayComp)
unique(bbd$Day) #"D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
bbd$Day = factor(bbd$Day, levels=c("D1","D3","D7","D10","D14","D21","D28","D36","D42"))
levels(bbd$Day) #"D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
(bbd_logfoldplot <- ggplot(data=bbd, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 5, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
bbd_logfoldplot


#BB and control, >1% abundance, at least 100 reads, select genera Log2fold plots
bbe <- read.csv('SRD129_FinalDiffAbundNasalGenus_BBControl_NoDNEG_CondensedList_GreaterThan1PercentAbundance.csv', header = TRUE, sep = ",")
head(bbe[,1:10])
colnames(bbe)
bbe$DayComp <- sub('_[A-Za-z]+', '\\2', bbe$comp)
unique(bbe$DayComp) #"D1_BBvsControl"  "D10_BBvsControl" "D3_BBvsControl"  "D7_BBvsControl"  "D14_BBvsControl" "D28_BBvsControl" "D42_BBvsControl" "D21_BBvsControl" "D36_BBvsControl"
bbe$Day <- sub('_[A-Za-z]+', '\\2', bbe$DayComp)
unique(bbe$Day) #"D1"  "D10" "D3"  "D7"  "D14" "D28" "D42" "D21" "D36"
bbe$Day = factor(bbe$Day, levels=c("D1","D3","D7","D10","D14","D21","D28","D36","D42"))
levels(bbe$Day) #"D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
(bbe_logfoldplot <- ggplot(data=bbe, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 5, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
bbe_logfoldplot <- bbe_logfoldplot + theme(strip.text = element_text(size= 13, face='italic'))
bbe_logfoldplot

#BB and control, >1% and <1% (respiratory organisms including Streptococcus, Actinobacillus, Acinetobacter, and Lactobacillus) abundance, at least 100 reads, select genera Log2fold plots
bbf <- read.csv('SRD129_FinalDiffAbundNasalGenus_BBControl_NoDNEG_CondensedList_GreaterThan1PercentAbundance_WithFS1RespiratoryOrganisms.csv', header = TRUE, sep = ",")
head(bbf[,1:10])
colnames(bbf)
bbf$DayComp <- sub('_[A-Za-z]+', '\\2', bbf$comp)
unique(bbf$DayComp) #"D1_BBvsControl"  "D10_BBvsControl" "D3_BBvsControl"  "D7_BBvsControl"  "D14_BBvsControl" "D28_BBvsControl" "D42_BBvsControl" "D21_BBvsControl" "D36_BBvsControl"
bbf$Day <- sub('_[A-Za-z]+', '\\2', bbf$DayComp)
unique(bbf$Day) #"D1"  "D10" "D3"  "D7"  "D14" "D28" "D42" "D21" "D36"
bbf$Day = factor(bbf$Day, levels=c("D1","D3","D7","D10","D14","D21","D28","D36","D42"))
levels(bbf$Day) #"D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
(bbf_logfoldplot <- ggplot(data=bbf, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 5, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
bbf_logfoldplot <- bbf_logfoldplot + theme(strip.text = element_text(size= 13, face='italic'))
bbf_logfoldplot

######################################################################################################

#BB between days Log2fold and Basemean plots
bb <- read.csv('SRD129_FinalDiffAbundNasalGenus_SignificantDaysWithinBBGroup.csv', header = TRUE, sep = ",")
head(bb[,1:10])
colnames(bb)
bb$DayComp <- sub('_[A-Za-z]+', '\\2', bb$comp)
unique(bb$DayComp) #"BB_D1vsD0"  "BB_D3vsD0"  "BB_D7vsD0"  "BB_D10vsD0" "BB_D14vsD0" "BB_D21vsD0" "BB_D28vsD0" "BB_D36vsD0" "BB_D42vsD0"
bb$Day <- sub('[A-Za-z]+_', '\\2', bb$DayComp)
unique(bb$Day) #"D1vsD0"  "D3vsD0"  "D7vsD0"  "D10vsD0" "D14vsD0" "D21vsD0" "D28vsD0" "D36vsD0" "D42vsD0"
bb$Day = factor(bb$Day, levels=c("D1vsD0",  "D3vsD0",  "D7vsD0",  "D10vsD0", "D14vsD0", "D21vsD0", "D28vsD0", "D36vsD0", "D42vsD0"))
levels(bb$Day)


(logfoldplot <- ggplot(data=bb, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 2)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
logfoldplot

(basemeanplot <- ggplot(data=bb, aes(x=Day, y=baseMean, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, scales = "free") + ylab('baseMean') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 2)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
basemeanplot


#BB between days select genera Log2fold plots
bbc <- read.csv('SRD129_FinalDiffAbundNasalSelectGenus_SignificantDaysWithinBBGroup.csv', header = TRUE, sep = ",")
head(bbc[,1:10])
colnames(bbc)
bbc$DayComp <- sub('_[A-Za-z]+', '\\2', bbc$comp)
unique(bbc$DayComp) #"BB_D1vsD0"  "BB_D3vsD0"  "BB_D7vsD0"  "BB_D10vsD0" "BB_D14vsD0" "BB_D21vsD0" "BB_D28vsD0" "BB_D36vsD0" "BB_D42vsD0"
bbc$Day <- sub('[A-Za-z]+_', '\\2', bbc$DayComp)
unique(bbc$Day) #"D1vsD0"  "D3vsD0"  "D7vsD0"  "D10vsD0" "D14vsD0" "D21vsD0" "D28vsD0" "D36vsD0" "D42vsD0"
bbc$Day = factor(bbc$Day, levels=c("D1vsD0",  "D3vsD0",  "D7vsD0",  "D10vsD0", "D14vsD0", "D21vsD0", "D28vsD0", "D36vsD0", "D42vsD0"))
levels(bbc$Day) #"D1vsD0"  "D3vsD0"  "D7vsD0"  "D10vsD0" "D14vsD0" "D21vsD0" "D28vsD0" "D36vsD0" "D42vsD0"
(bbc_logfoldplot <- ggplot(data=bbc, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 7, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
bbc_logfoldplot


#BB between days even smaller list of genera Log2fold plots
bbb <- read.csv('SRD129_FinalDiffAbundNasalSmallerListGenus_SignificantDaysWithinBBGroup.csv', header = TRUE, sep = ",")
head(bbb[,1:10])
colnames(bbb)
bbb$DayComp <- sub('_[A-Za-z]+', '\\2', bbb$comp)
unique(bbb$DayComp) #"BB_D1vsD0"  "BB_D3vsD0"  "BB_D7vsD0"  "BB_D10vsD0" "BB_D14vsD0" "BB_D21vsD0" "BB_D28vsD0" "BB_D36vsD0" "BB_D42vsD0"
bbb$Day <- sub('[A-Za-z]+_', '\\2', bbb$DayComp)
unique(bbb$Day) #"D1vsD0"  "D3vsD0"  "D7vsD0"  "D10vsD0" "D14vsD0" "D21vsD0" "D28vsD0" "D36vsD0" "D42vsD0"
bbb$Day = factor(bbb$Day, levels=c("D1vsD0",  "D3vsD0",  "D7vsD0",  "D10vsD0", "D14vsD0", "D21vsD0", "D28vsD0", "D36vsD0", "D42vsD0"))
levels(bbb$Day) #"D1vsD0"  "D3vsD0"  "D7vsD0"  "D10vsD0" "D14vsD0" "D21vsD0" "D28vsD0" "D36vsD0" "D42vsD0"
(bbb_logfoldplot <- ggplot(data=bbb, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 7, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
bbb_logfoldplot





######################################################################################################

unique(final.nonsigtab$OTU) #number of unique OTUs in all groups combined
#  [1] "Otu0020" "Otu0087" "Otu0115" "Otu0490" "Otu0025" "Otu0015" "Otu0016" "Otu0022" "Otu0024" "Otu0026" "Otu0027" "Otu0029" "Otu0033"
# [14] "Otu0034" "Otu0041" "Otu0049" "Otu0058" "Otu0060" "Otu0063" "Otu0066" "Otu0070" "Otu0074" "Otu0079" "Otu0082" "Otu0083" "Otu0088"
# [27] "Otu0089" "Otu0097" "Otu0104" "Otu0114" "Otu0122" "Otu0123" "Otu0126" "Otu0128" "Otu0131" "Otu0134" "Otu0151" "Otu0156" "Otu0157"
# [40] "Otu0162" "Otu0169" "Otu0174" "Otu0176" "Otu0177" "Otu0183" "Otu0184" "Otu0209" "Otu0232" "Otu0253" "Otu0256" "Otu0263" "Otu0268"
# [53] "Otu0278" "Otu0293" "Otu0294" "Otu0310" "Otu0320" "Otu0326" "Otu0348" "Otu0359" "Otu0360" "Otu0376" "Otu0377" "Otu0392" "Otu0402"
# [66] "Otu0407" "Otu0409" "Otu0420" "Otu0427" "Otu0435" "Otu0470" "Otu0480" "Otu0520" "Otu0635" "Otu0674" "Otu0733" "Otu0975" "Otu1003"
# [79] "Otu0001" "Otu0006" "Otu0017" "Otu0021" "Otu0028" "Otu0047" "Otu0051" "Otu0053" "Otu0054" "Otu0059" "Otu0072" "Otu0075" "Otu0080"
# [92] "Otu0086" "Otu0100" "Otu0102" "Otu0110" "Otu0127" "Otu0135" "Otu0150" "Otu0167" "Otu0192" "Otu0194" "Otu0196" "Otu0200" "Otu0203"
# [105] "Otu0221" "Otu0244" "Otu0255" "Otu0260" "Otu0262" "Otu0265" "Otu0292" "Otu0342" "Otu0384" "Otu0386" "Otu0406" "Otu0438" "Otu0472"
# [118] "Otu0473" "Otu0485" "Otu0565" "Otu0573" "Otu0672" "Otu0688" "Otu0745" "Otu0949" "Otu0952" "Otu0030" "Otu0056" "Otu0068" "Otu0076"
# [131] "Otu0091" "Otu0149" "Otu0154" "Otu0297" "Otu0298" "Otu0313" "Otu0340" "Otu0351" "Otu0363" "Otu0404" "Otu0525" "Otu0529" "Otu0701"
# [144] "Otu0752" "Otu0941" "Otu0004" "Otu0013" "Otu0018" "Otu0037" "Otu0057" "Otu0073" "Otu0108" "Otu0129" "Otu0216" "Otu0229" "Otu0236"
# [157] "Otu0246" "Otu0259" "Otu0283" "Otu0486" "Otu0032" "Otu0035" "Otu0185" "Otu0234" "Otu0299" "Otu0322" "Otu0398" "Otu0537" "Otu0009"
# [170] "Otu0061" "Otu0065" "Otu0084" "Otu0269" "Otu0457" "Otu0423" "Otu0046" "Otu0067" "Otu0188" "Otu0218" "Otu0393" "Otu0005" "Otu0011"
# [183] "Otu0012" "Otu0350" "Otu0602" "Otu0038" "Otu0193" "Otu0343" "Otu0371" "Otu0580" "Otu0677" "Otu0036" "Otu0224" "Otu0724" "Otu0002"
# [196] "Otu0116" "Otu0170" "Otu0304" "Otu0362" "Otu0477" "Otu0467" "Otu0266" "Otu0498" "Otu0144" "Otu0284" "Otu0646" "Otu0605" "Otu0714"
# [209] "Otu0287" "Otu0288" "Otu0330" "Otu0347" "Otu0432" "Otu0487" "Otu0558" "Otu0649" "Otu0766" "Otu0878" "Otu0338" "Otu0475" "Otu0719"
# [222] "Otu0793" "Otu0832" "Otu0462" "Otu0837" "Otu0315" "Otu0433" "Otu0452" "Otu0514" "Otu0535" "Otu0665" "Otu0308" "Otu0629" "Otu0173"
# [235] "Otu0385" "Otu0499" "Otu0327" "Otu0775" "Otu0399" "Otu0551" "Otu0631" "Otu0698" "Otu0702" "Otu0917" "Otu0413" "Otu0458" "Otu0508"
# [248] "Otu0604" "Otu0621" "Otu0715" "Otu0739" "Otu0771" "Otu0777" "Otu0859" "Otu1260" "Otu0934" "Otu0765" "Otu0922" "Otu1036" "Otu1090"
# [261] "Otu0093"

#261 unique OTUs from non-significant days

#######################################################################################################

#Generate a plot of color names which R knows about: http://www.sthda.com/english/wiki/colors-in-r
showCols <- function(cl=colors(), bg = "grey",
                     cex = 0.75, rot = 30) {
  m <- ceiling(sqrt(n <-length(cl)))
  length(cl) <- m*m; cm <- matrix(cl, m)
  require("grid")
  grid.newpage(); vp <- viewport(w = .92, h = .92)
  grid.rect(gp=gpar(fill=bg))
  grid.text(cm, x = col(cm)/m, y = rev(row(cm))/m, rot = rot,
            vp=vp, gp=gpar(cex = cex, col = cm))
}


showCols(cl = colors(), bg= "gray33", rot=30, cex=0.75)
