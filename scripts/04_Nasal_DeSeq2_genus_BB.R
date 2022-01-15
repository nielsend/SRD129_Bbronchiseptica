#####################################################################################################
#SRD129 16S - Differential Abundance of Genera in Nasal Microbiota using DESeq2
#By Mou, KT

#NOTES: 
#This code analyzes differential abundance of OTUs in nasal samples between BB and control groups
#using DESeq2 at the genus level

#Clear workspace and load necessary packages
rm(list=ls())

#Files needed:
#Mothur shared file: SRD129BB.outsingletons.abund.opti_mcc.shared
#Mothur constaxonomy file: SRD129BB.outsingletons.abund.opti_mcc.0.03.cons.taxonomy
#Metadata: SRD129BBmetadata.csv

#Load library packages
library(DESeq2)
library(phyloseq)
library(ggplot2)
library(tidyr)
library("wesanderson")
library(plotly)
library(gapminder)
library("ggsci")

#Set plots to have gray background with white gridlines
theme_set(theme_gray())

########################################################################################################

####### PREPARING OBJECTS FOR DESEQ2 ANALYSIS ########

#Load files
otu2 <- import_mothur(mothur_shared_file = './data/SRD129BB.outsingletons.abund.opti_mcc.shared') #use unrarified data
taxo2 <- import_mothur(mothur_constaxonomy_file = './data/SRD129BB.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
taxo2
meta2 <- read.table(file = './data/SRD129BBmetadata.csv', sep = ',', header = TRUE)

#Organize meta file
rownames(meta2) <- meta2$Sample
meta2 <- meta2[,-1] #remove Sample column
meta2$Set <- paste(meta2$Day, meta2$Treatment, sep = '_')

#Make phyloseq object SRD129 (combine taxonomy, OTU, and metadata)
phy_meta2 <- sample_data(meta2) 
SRD129 <- phyloseq(otu2, taxo2)
SRD129 <- merge_phyloseq(SRD129, phy_meta2)   # combines the metadata with this phyloseq object
colnames(tax_table(SRD129))
colnames(tax_table(SRD129)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
SRD129
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1947 taxa and 190 samples ]
# sample_data() Sample Data:       [ 190 samples by 5 sample variables ]
# tax_table()   Taxonomy Table:    [ 1947 taxa by 6 taxonomic ranks ]

#Prune
SRD129 <- prune_samples(sample_sums(SRD129) > 2000, SRD129)  # This removes samples that have fewer than 2000 sequences associated with them.
SRD129 <- prune_taxa(taxa_sums(SRD129) > 10, SRD129)        # removes OTUs that occur less than 10 times globally
tax_table(SRD129) [1:5, 1:6] #see what's in tax_table first 5 rows, first 6 columns

# If you want to group OTUs uncomment the tax_glom() line and select your desired taxrank
# right now all these analysis are done at the OTU level.

SRD129.genus <- tax_glom(SRD129, taxrank = "Genus")
# This method merges species that have the same taxonomy at a certain taxanomic rank. 
# Its approach is analogous to tip_glom, but uses categorical data instead of a tree. 

#Subset nasal samples to count number of sequences for each OTU
SRD129.nw <- psmelt(SRD129)
write.csv(SRD129.nw, file="SRD129BB_GenusSampleData.csv")

SRD129.nw.otu <- otu_table(SRD129)
write.csv(SRD129.nw.otu, file="SRD129_DESeqOTUtable.csv")

########################### Day 0 #########################################################

sample_data(SRD129.genus)

SRD129.D0 <- subset_samples(SRD129.genus, Day == 'D0')
sample_sums(SRD129.D0)
colnames(otu_table(SRD129.D0)) #check on all the sample names
SRD129.D0 <- prune_taxa(taxa_sums(SRD129.D0) > 1, SRD129.D0)
#if taxa_sums is >1, then it will print that out in SRD129.D0 object and not include anything with <1.
rowSums(SRD129.D0@otu_table)
SRD129.D0.De <- phyloseq_to_deseq2(SRD129.D0, ~ Set)
# ~Set: whatever you want to group data by, whatever column you used to designate ellipses with

#DESeq calculation: Differential expression analysis 
#based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
SRD129.D0.De <- DESeq(SRD129.D0.De, test = "Wald", fitType = "parametric")

######### Day 0 BB vs Control ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D0_BB")
#BB = 10
sum(meta2$Set == "D0_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D0.bc = results(SRD129.D0.De, contrast=c("Set","D0_BB","D0_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#used pAdjust method = BH because that was what the tutorial used
res.D0.bc
sigtab.D0.bc = res.D0.bc[which(res.D0.bc$padj < .05), ]
sigtab.D0.bc = cbind(as(sigtab.D0.bc, "data.frame"), as(tax_table(SRD129.D0)[rownames(sigtab.D0.bc), ], "matrix"))
#cbind combines all columns together, regardless of rownames (match by rownames, use merge function)
format(sigtab.D0.bc$padj, scientific = TRUE)
sigtab.D0.bc$newp <- format(round(sigtab.D0.bc$padj, digits = 3), scientific = TRUE)
sigtab.D0.bc$Treatment <- ifelse(sigtab.D0.bc$log2FoldChange >=0, "BB", "Control")
#ifelse: Assigning "BB" = yes, "Control" = no; it's important to make sure you have the
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
sigtab.D0.bc$comp <- 'D0_BBvsControl'

#Create final significant comparisons table
final.nonsigtab <- sigtab.D0.bc

############## Day 1 #########################

sample_data(SRD129.genus)
SRD129.D1 <- subset_samples(SRD129.genus, Day == 'D1')
sample_sums(SRD129.D1)
colnames(otu_table(SRD129.D1)) #check on all the sample names
SRD129.D1 <- prune_taxa(taxa_sums(SRD129.D1) > 1, SRD129.D1)
rowSums(SRD129.D1@otu_table)
SRD129.D1.De <- phyloseq_to_deseq2(SRD129.D1, ~ Set)
SRD129.D1.De <- DESeq(SRD129.D1.De, test = "Wald", fitType = "parametric")

######### Day 1 BB vs Control ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D1_BB")
#BB = 9
sum(meta2$Set == "D1_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D1.bc = results(SRD129.D1.De, contrast=c("Set","D1_BB","D1_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D1.bc
sigtab.D1.bc = res.D1.bc[which(res.D1.bc$padj < .05), ]
sigtab.D1.bc = cbind(as(sigtab.D1.bc, "data.frame"), as(tax_table(SRD129.D1)[rownames(sigtab.D1.bc), ], "matrix"))
format(sigtab.D1.bc$padj, scientific = TRUE)
sigtab.D1.bc$newp <- format(round(sigtab.D1.bc$padj, digits = 3), scientific = TRUE)
sigtab.D1.bc$Treatment <- ifelse(sigtab.D1.bc$log2FoldChange >=0, "BB", "Control")
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
sigtab.D1.bc$comp <- 'D1_BBvsControl'

#Add to final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D1.bc)

############## Day 3 #########################

sample_data(SRD129.genus)
SRD129.D3 <- subset_samples(SRD129.genus, Day == 'D3')
sample_sums(SRD129.D3)
colnames(otu_table(SRD129.D3)) #check on all the sample names
SRD129.D3 <- prune_taxa(taxa_sums(SRD129.D3) > 1, SRD129.D3)
rowSums(SRD129.D3@otu_table)
SRD129.D3.De <- phyloseq_to_deseq2(SRD129.D3, ~ Set)
SRD129.D3.De <- DESeq(SRD129.D3.De, test = "Wald", fitType = "parametric")

######### Day 3 BB vs Control ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D3_BB")
#BB = 9
sum(meta2$Set == "D3_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D3.bc = results(SRD129.D3.De, contrast=c("Set","D3_BB","D3_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D3.bc
sigtab.D3.bc = res.D3.bc[which(res.D3.bc$padj < .05), ]
sigtab.D3.bc = cbind(as(sigtab.D3.bc, "data.frame"), as(tax_table(SRD129.D3)[rownames(sigtab.D3.bc), ], "matrix"))
format(sigtab.D3.bc$padj, scientific = TRUE)
sigtab.D3.bc$newp <- format(round(sigtab.D3.bc$padj, digits = 3), scientific = TRUE)
sigtab.D3.bc$Treatment <- ifelse(sigtab.D3.bc$log2FoldChange >=0, "BB", "Control")
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
sigtab.D3.bc$comp <- 'D3_BBvsControl'

#Add to final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D3.bc)

############## Day 7 #########################

sample_data(SRD129.genus)
SRD129.D7 <- subset_samples(SRD129.genus, Day == 'D7')
sample_sums(SRD129.D7)
colnames(otu_table(SRD129.D7)) #check on all the sample names
SRD129.D7 <- prune_taxa(taxa_sums(SRD129.D7) > 1, SRD129.D7)
rowSums(SRD129.D7@otu_table)
SRD129.D7.De <- phyloseq_to_deseq2(SRD129.D7, ~ Set)
SRD129.D7.De <- DESeq(SRD129.D7.De, test = "Wald", fitType = "parametric")

######### Day 7 BB vs Control ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D7_BB")
#BB = 9
sum(meta2$Set == "D7_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D7.bc = results(SRD129.D7.De, contrast=c("Set","D7_BB","D7_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D7.bc
sigtab.D7.bc = res.D7.bc[which(res.D7.bc$padj < .05), ]
sigtab.D7.bc = cbind(as(sigtab.D7.bc, "data.frame"), as(tax_table(SRD129.D7)[rownames(sigtab.D7.bc), ], "matrix"))
format(sigtab.D7.bc$padj, scientific = TRUE)
sigtab.D7.bc$newp <- format(round(sigtab.D7.bc$padj, digits = 3), scientific = TRUE)
sigtab.D7.bc$Treatment <- ifelse(sigtab.D7.bc$log2FoldChange >=0, "BB", "Control")
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
sigtab.D7.bc$comp <- 'D7_BBvsControl'

#Add to final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D7.bc)

############## Day 10 #########################

sample_data(SRD129.genus)
SRD129.D10 <- subset_samples(SRD129.genus, Day == 'D10')
sample_sums(SRD129.D10)
colnames(otu_table(SRD129.D10)) #check on all the sample names
SRD129.D10 <- prune_taxa(taxa_sums(SRD129.D10) > 1, SRD129.D10)
rowSums(SRD129.D10@otu_table)
SRD129.D10.De <- phyloseq_to_deseq2(SRD129.D10, ~ Set)
SRD129.D10.De <- DESeq(SRD129.D10.De, test = "Wald", fitType = "parametric")

######### Day 10 BB vs Control ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D10_BB")
#BB = 10
sum(meta2$Set == "D10_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D10.bc = results(SRD129.D10.De, contrast=c("Set","D10_BB","D10_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D10.bc
sigtab.D10.bc = res.D10.bc[which(res.D10.bc$padj < .05), ]
sigtab.D10.bc = cbind(as(sigtab.D10.bc, "data.frame"), as(tax_table(SRD129.D10)[rownames(sigtab.D10.bc), ], "matrix"))
format(sigtab.D10.bc$padj, scientific = TRUE)
sigtab.D10.bc$newp <- format(round(sigtab.D10.bc$padj, digits = 3), scientific = TRUE)
sigtab.D10.bc$Treatment <- ifelse(sigtab.D10.bc$log2FoldChange >=0, "BB", "Control")
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
sigtab.D10.bc$comp <- 'D10_BBvsControl'

#Add to final significant comparisons table of days
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D10.bc)

############## Day 14 #########################

sample_data(SRD129.genus)
SRD129.D14 <- subset_samples(SRD129.genus, Day == 'D14')
sample_sums(SRD129.D14)
colnames(otu_table(SRD129.D14)) #check on all the sample names
SRD129.D14 <- prune_taxa(taxa_sums(SRD129.D14) > 1, SRD129.D14)
rowSums(SRD129.D14@otu_table)
SRD129.D14.De <- phyloseq_to_deseq2(SRD129.D14, ~ Set)
SRD129.D14.De <- DESeq(SRD129.D14.De, test = "Wald", fitType = "parametric")

######### Day 14 BB vs Control ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D14_BB")
#BB = 10
sum(meta2$Set == "D14_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D14.bc = results(SRD129.D14.De, contrast=c("Set","D14_BB","D14_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D14.bc
sigtab.D14.bc = res.D14.bc[which(res.D14.bc$padj < .05), ]
sigtab.D14.bc = cbind(as(sigtab.D14.bc, "data.frame"), as(tax_table(SRD129.D14)[rownames(sigtab.D14.bc), ], "matrix"))
format(sigtab.D14.bc$padj, scientific = TRUE)
sigtab.D14.bc$newp <- format(round(sigtab.D14.bc$padj, digits = 3), scientific = TRUE)
sigtab.D14.bc$Treatment <- ifelse(sigtab.D14.bc$log2FoldChange >=0, "BB", "Control")
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
sigtab.D14.bc$comp <- 'D14_BBvsControl'

#Add to final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D14.bc)

############## Day 21 #########################

sample_data(SRD129.genus)
SRD129.D21 <- subset_samples(SRD129.genus, Day == 'D21')
sample_sums(SRD129.D21)
colnames(otu_table(SRD129.D21)) #check on all the sample names
SRD129.D21 <- prune_taxa(taxa_sums(SRD129.D21) > 1, SRD129.D21)
rowSums(SRD129.D21@otu_table)
SRD129.D21.De <- phyloseq_to_deseq2(SRD129.D21, ~ Set)
SRD129.D21.De <- DESeq(SRD129.D21.De, test = "Wald", fitType = "parametric")

######### Day 21 BB vs Control ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D21_BB")
#BB = 10
sum(meta2$Set == "D21_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D21.bc = results(SRD129.D21.De, contrast=c("Set","D21_BB","D21_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D21.bc
sigtab.D21.bc = res.D21.bc[which(res.D21.bc$padj < .05), ]
sigtab.D21.bc = cbind(as(sigtab.D21.bc, "data.frame"), as(tax_table(SRD129.D21)[rownames(sigtab.D21.bc), ], "matrix"))
format(sigtab.D21.bc$padj, scientific = TRUE)
sigtab.D21.bc$newp <- format(round(sigtab.D21.bc$padj, digits = 3), scientific = TRUE)
sigtab.D21.bc$Treatment <- ifelse(sigtab.D21.bc$log2FoldChange >=0, "BB", "Control")
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
sigtab.D21.bc$comp <- 'D21_BBvsControl'

#Add to final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D21.bc)

############## Day 28 #########################

sample_data(SRD129.genus)
SRD129.D28 <- subset_samples(SRD129.genus, Day == 'D28')
sample_sums(SRD129.D28)
colnames(otu_table(SRD129.D28)) #check on all the sample names
SRD129.D28 <- prune_taxa(taxa_sums(SRD129.D28) > 1, SRD129.D28)
rowSums(SRD129.D28@otu_table)
SRD129.D28.De <- phyloseq_to_deseq2(SRD129.D28, ~ Set)
SRD129.D28.De <- DESeq(SRD129.D28.De, test = "Wald", fitType = "parametric")

######### Day 28 BB vs Control ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D28_BB")
#BB = 10
sum(meta2$Set == "D28_Control")
#Control = 4

#Extract results from a DESeq analysis, organize table
res.D28.bc = results(SRD129.D28.De, contrast=c("Set","D28_BB","D28_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D28.bc
sigtab.D28.bc = res.D28.bc[which(res.D28.bc$padj < .05), ]
sigtab.D28.bc = cbind(as(sigtab.D28.bc, "data.frame"), as(tax_table(SRD129.D28)[rownames(sigtab.D28.bc), ], "matrix"))
format(sigtab.D28.bc$padj, scientific = TRUE)
sigtab.D28.bc$newp <- format(round(sigtab.D28.bc$padj, digits = 3), scientific = TRUE)
sigtab.D28.bc$Treatment <- ifelse(sigtab.D28.bc$log2FoldChange >=0, "BB", "Control")
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
sigtab.D28.bc$comp <- 'D28_BBvsControl'

#Add to final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D28.bc)

############## Day 36 #########################

sample_data(SRD129.genus)
SRD129.D36 <- subset_samples(SRD129.genus, Day == 'D36')
sample_sums(SRD129.D36)
colnames(otu_table(SRD129.D36)) #check on all the sample names
SRD129.D36 <- prune_taxa(taxa_sums(SRD129.D36) > 1, SRD129.D36)
rowSums(SRD129.D36@otu_table)
SRD129.D36.De <- phyloseq_to_deseq2(SRD129.D36, ~ Set)
SRD129.D36.De <- DESeq(SRD129.D36.De, test = "Wald", fitType = "parametric")

######### Day 36 BB vs Control ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D36_BB")
#BB = 10
sum(meta2$Set == "D36_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D36.bc = results(SRD129.D36.De, contrast=c("Set","D36_BB","D36_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D36.bc
sigtab.D36.bc = res.D36.bc[which(res.D36.bc$padj < .05), ]
sigtab.D36.bc = cbind(as(sigtab.D36.bc, "data.frame"), as(tax_table(SRD129.D36)[rownames(sigtab.D36.bc), ], "matrix"))
format(sigtab.D36.bc$padj, scientific = TRUE)
sigtab.D36.bc$newp <- format(round(sigtab.D36.bc$padj, digits = 3), scientific = TRUE)
sigtab.D36.bc$Treatment <- ifelse(sigtab.D36.bc$log2FoldChange >=0, "BB", "Control")
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
sigtab.D36.bc$comp <- 'D36_BBvsControl'

#Create final significant comparisons table
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D36.bc)

############## Day 42 #########################

sample_data(SRD129.genus)
SRD129.D42 <- subset_samples(SRD129.genus, Day == 'D42')
sample_sums(SRD129.D42)
colnames(otu_table(SRD129.D42)) #check on all the sample names
SRD129.D42 <- prune_taxa(taxa_sums(SRD129.D42) > 1, SRD129.D42)
rowSums(SRD129.D42@otu_table)
SRD129.D42.De <- phyloseq_to_deseq2(SRD129.D42, ~ Set)
SRD129.D42.De <- DESeq(SRD129.D42.De, test = "Wald", fitType = "parametric")

######### Day 42 BB vs Control ###################

#Number of pigs per group (using meta2 dataframe): 
sum(meta2$Set == "D42_BB")
#BB = 10
sum(meta2$Set == "D42_Control")
#Control = 10

#Extract results from a DESeq analysis, organize table
res.D42.bc = results(SRD129.D42.De, contrast=c("Set","D42_BB","D42_Control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D42.bc
sigtab.D42.bc = res.D42.bc[which(res.D42.bc$padj < .05), ]
sigtab.D42.bc = cbind(as(sigtab.D42.bc, "data.frame"), as(tax_table(SRD129.D42)[rownames(sigtab.D42.bc), ], "matrix"))
format(sigtab.D42.bc$padj, scientific = TRUE)
sigtab.D42.bc$newp <- format(round(sigtab.D42.bc$padj, digits = 3), scientific = TRUE)
sigtab.D42.bc$Treatment <- ifelse(sigtab.D42.bc$log2FoldChange >=0, "BB", "Control")
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
sigtab.D42.bc$comp <- 'D42_BBvsControl'

#Create final significant comparisons table of days with non-significant beta diversity changes
final.nonsigtab <- rbind(final.nonsigtab,sigtab.D42.bc)

#write csv and save within ./results
write.csv(final.nonsigtab, file= "SRD129_BBControl_FinalDiffAbundNasalGenus_SignificantDays.csv")

#Export nasalgenplot from 05_Genus_BB.R as a figure and cross-reference with SRD129_BBControl_FinalDiffAbundNasalGenus_SignificantDays.csv
#to narrow down list of DESeq2 significant genera with abundance >1%, and are not differentially abundant on day 0.
#Save a copy of SRD129_BBControl_FinalDiffAbundNasalGenus_SignificantDays.csv with the narrowed list of genera and
#save as SRD129_BBControl_Nasal_GenusAbundanceMatchDESeq2List.csv

#Plot new list of DESeq2 significant genera
nasalgen2 <- read.csv("./data/SRD129_BBControl_Nasal_GenusAbundanceMatchDESeq2List.csv")
nasalgen2$Day <- sub('_[A-Za-z]+', '\\2', nasalgen2$comp) #create new column "Day" by extracting from column "comp"
unique(nasalgen2$Day) #"D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
nasalgen2$Day = factor(nasalgen2$Day, levels=c("D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
levels(nasalgen2$Day) #"D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
unique(nasalgen2$Treatment) #BB Control
(nasalgen2_logfoldplot <- ggplot(data=nasalgen2, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, ncol = 5, scales = "free") + ylab('log2-fold change') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
(nasalgen2_logfoldplot <- nasalgen2_logfoldplot + 
    theme(strip.text = element_text(size= 13, face='italic'), 
          legend.text = element_text(size=13),
          legend.title = element_text(size=13)))

#Save 'nasalgen2_logfoldplot' as a .tiff for publication, 500dpi
ggsave("Figure_5.tiff", plot=nasalgen2_logfoldplot, width = 16, height = 6, dpi = 500, units =c("in"))