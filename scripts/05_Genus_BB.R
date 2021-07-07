############################################################
#SRD129 16S - Total Genera
#By Mou KT

#NOTES: 
#This code determines the total number of genera found in each treatment group per day and plots as bar graphs
#It also assesses total number of unique families in each day+treatment group 

#Clear workspace and load necessary packages
rm(list=ls())

#Files needed:
#Mothur subsample.shared file: SRD129BB.outsingletons.abund.opti_mcc.0.03.subsample.shared
#Mothur constaxonomy file: SRD129BB.outsingletons.abund.opti_mcc.0.03.cons.taxonomy
#Metadata: SRD129BBmetadata.csv

#Load libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(cowplot)
library("ggsci")

############################################################

otu <- import_mothur(mothur_shared_file = './data/SRD129BB.outsingletons.abund.opti_mcc.0.03.subsample.shared')
sample_sums(otu)
taxo <- import_mothur(mothur_constaxonomy_file = './data/SRD129BB.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
shared <- read.table('./data/SRD129BB.outsingletons.abund.opti_mcc.0.03.subsample.shared', header = TRUE)
meta <- read.table('./data/SRD129BBmetadata.csv', header = TRUE, sep = ",")
head(meta)

colnames(meta)[1] <- 'group' #rename as group temporarily. Will use to set as rownames later and remove the
#"group" column
#meta$Day<- gsub("D", "", meta$Day)
meta$group <- as.character(meta$group)
head(meta)

phy_meta <- sample_data(meta) 
rownames(phy_meta) <- phy_meta$group #Set group names as rownames
head(phy_meta)
phy_meta <- phy_meta[,-1] #remove group column
head(phy_meta)
sample_sums(otu)

SRD129 <- phyloseq(otu, taxo) #build phyloseq-class objects from their components
SRD129 <- merge_phyloseq(SRD129, phy_meta)  # combines the metadata with this phyloseq object
colnames(tax_table(SRD129)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
sample_sums(SRD129) #sum of all OTUs for each sample, almost all samples have 2000 sequences
SRD129 <- prune_taxa(taxa_sums(SRD129) > 2, SRD129)  # removes OTUs that occur less than 2 times globally
SRD129.genus <- tax_glom(SRD129, 'Genus')
phyla_tab <- as.data.frame(t(SRD129.genus@otu_table)) #transpose all
head(phyla_tab)
SRD129.genus@tax_table[,6]
colnames(phyla_tab) <- SRD129.genus@tax_table[,6] #replace column names in phyla_tab from Otuxxxx with Genus types
#5 refers to second column of tax_table(SRD129) which is 'Genus'
phyla_tab2 <- phyla_tab/rowSums(phyla_tab) #proportion of specific phyla per phyla column in phyla_tab
head(phyla_tab2)
phyla_tab2$group <- rownames(phyla_tab2) #create new column called "group" containing rownames
head(phyla_tab2)
fobar <- merge(meta, phyla_tab2, by = 'group') #merge meta with phyla_tab2 by "group"
head(fobar)

fobar.gather <- fobar %>% gather(Genus, value, -(group:Treatment))  #converts to long form dataframe, this is handy for using ggplot faceting functions, check out tidyverse tutorials
#create new columns Genus, value; add columns group to All before Genus and value
head(fobar.gather)

#modified ggplot
ggplot(data=fobar.gather, aes(x=Treatment, y=value, fill=Genus)) +
  geom_bar(stat = 'identity') +
  facet_grid(Treatment~Day) + ylab('Percent of total genus in BB and Control groups')

#check if there is an extra "group" column. If so, run the next command and modify which column to remove
which(colnames(phyla_tab2)=="group") #result says column 288
phyla_tab3 <- phyla_tab2[,-288] ##drop the 288th column because it's an extra group (sample name) column
phyla_tab4 <- phyla_tab3[,colSums(phyla_tab3)>0.1] #keep the columns that have greater than 0.1 value
phyla_tab4$group <- rownames(phyla_tab4) #rename rownames as group
fobar2 <- merge(meta, phyla_tab4, by = 'group')
head(fobar2)
fobar2.gather <- fobar2 %>% gather(Genus, value, -(group:Treatment))  #converts to long form dataframe, this is handy for using ggplot faceting functions, check out tidyverse tutorials
head(fobar2.gather)
fobar2.gather$All <- paste(fobar2.gather$Day, fobar2.gather$Treatment, sep = '_')
#Create "All" column with Day, Treatment
#Use this to make % total community plot
head(fobar2.gather$All)

#Ways to summarize unique genera
unique(fobar2.gather$Genus) #Total number of unique genera: 95
#OR
fobar2.gather %>% summarise_each(funs(n_distinct)) #95 total unique genus

#Calculate percent abundance of each genus per "All" type
fobar2.gather <- fobar2.gather %>% group_by(All) %>% mutate(value2=(value/(length(All)/95))*100)
#95 refers to number of Genus (counted based on summarise_each function that was ran on fobar2.gather
#to get total number unique items under "Genus" column)

#Graph percent abundance of each genus
theme_update(plot.title = element_text(hjust = 0.5)) #center title
ggplot(data=fobar2.gather, aes(x=Treatment, y=value2, fill=Genus)) +
  geom_bar(stat = 'identity') +
  facet_grid(Treatment~Day) + ylab('Percent of total community') +
  ggtitle("Total percent genus members in nasal wash \namong treatment groups")

write.csv(fobar2.gather, file="SRD129_BB_NasalGenusPercentAbundance.csv")

#############NASAL#############   CONTINUE WORKING FROM HERE!
#Subset each "All" group from fobar2.gather and save as individual csv files.
#For example:
D28_Control.genus <- subset(fobar2.gather, All=="D28_Control") #EDIT THIS
write.csv(D28_Control.genus, file="D28_Control.genus.csv")  #EDIT THIS
D28_BB.genus <- subset(fobar2.gather, All=="D28_BB")
write.csv(D28_BB.genus, file="D28_BB.genus.csv")

#Calculate the total % percent abundance of each genera on each sample (I used JMP to do this) 
#In JMP: copy and paste values from Genus and value columns
#Tables > Summary > input columns for Statistics (% of Total of value column) and Group (Genus column) > OK
#Copy results to a spreadsheet editor (see D0_Control.genus.xlsx for an example)
#Since we are only interested in genera that are above 1% abundance, 
#calculate total percentage of all other genera that are less than 1% in the spreadsheet and label as "Other".
#Create a new spreadsheet and label as "Nasal_genus.csv"
#Create the following columns: Day, Treatment, Genus, Percent Abundance
#Copy the list of genera and their percent abundances from each of the individual Excel files to "Nasal genus.csv" spreadsheet.
#Fill in the other columns manually (Day, Treatment). 
#You should have a file similar to SRD129_BBControl_GenusPercentAbundance.csv. Continue to the next step.

#Jan. 30, 2019: edit this csv file to include updated D1 BB percentage info
nasalgen <- read.table('SRD129_Nasal_GenusPercentAbundance.csv', header = TRUE, sep = ",")
head(nasalgen)



#####Nasal genus with other category above 1% abundance######
#Nasalgen.2 (BB and control only, more than 1% genera)
unique(nasalgen.2$Treatment) #Control BB
unique(nasalgen.2$Day) #D0  D1  D10 D14 D21 D28 D3  D36 D42 D7
nasalgen.2$Day = factor(nasalgen.2$Day, levels = c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
nasalgen.2$More.than.2=as.character(nasalgen.2$Genus)
str(nasalgen.2$More.than.2) #Compactly display the internal structure of an R object,
nasalgen.2$More.than.2[nasalgen.2$Percent.abundance<1]<-"Other"
write.csv(nasalgen.2, file = "SRD129_Nasal_BBControlNoDNEGGenusPercentAbundanceAbove1percent.csv")
nasalgen.2b <- read.table('SRD129_Nasal_BBControlNoDNEGGenusPercentAbundanceAbove1percentAddTo100.csv', header = TRUE, sep = ",")
head(nasalgen.2b)
levels(nasalgen.2b$Day) #"D0"  "D1"  "D10" "D14" "D21" "D28" "D3"  "D36" "D42" "D7" 
nasalgen.2b$Day = factor(nasalgen.2b$Day, levels = c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
levels(nasalgen.2b$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
write.csv(nasalgen.2b, file= "SRD129_BBControl_NasalGen.csv")

#Nasalgen.2 abundance plot for each genus, BB and control only, no DNEG12/6, more than 1% genera -- USE THIS
(nasalgen.2bplot2 <- ggplot(data=nasalgen.2b, aes(x=Day, y=Percent.abundance, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~More.than.2, scales = "free") + ylab('Relative abundance (%) at nasal site') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 2)) +
    theme(axis.line = element_line()) +
    scale_fill_igv(name = "Genus") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
nasalgen.2bplot2


#### Subset each "All" group from fobar2.gather ####
fobar2.gather$All
#Day 0 Control nasal total Genus
D0_Control_N.genus <- filter(fobar2.gather, grepl(paste("D0_Control"), All))
D0_Control_N.genus
D0_Control_N.genus$Indicator<-"Zero"
D0_Control_N.genus$Indicator[D0_Control_N.genus$value>0]<-D0_Control_N.genus$Genus[D0_Control_N.genus$value>0]
unique(D0_Control_N.genus$Indicator) #89 for indicator column
write.csv(D0_Control_N.genus, file = "D0_Control_N.genus.csv")

#Day 0 BB nasal total Genus
D0_BB_N.genus <- filter(fobar2.gather, grepl(paste("D0_BB"), All))
D0_BB_N.genus
D0_BB_N.genus$Indicator<-"Zero"
D0_BB_N.genus$Indicator[D0_BB_N.genus$value>0]<-D0_BB_N.genus$Genus[D0_BB_N.genus$value>0]
unique(D0_BB_N.genus$Indicator) #86 for indicator column
write.csv(D0_BB_N.genus, file = "D0_BB_N.genus.csv")

#Day 1 Control nasal total Genus
D1_Control_N.genus <- filter(fobar2.gather, grepl(paste("D1_Control"), All))
D1_Control_N.genus
D1_Control_N.genus$Indicator<-"Zero"
D1_Control_N.genus$Indicator[D1_Control_N.genus$value>0]<-D1_Control_N.genus$Genus[D1_Control_N.genus$value>0]
unique(D1_Control_N.genus$Indicator) #89 for indicator column
write.csv(D1_Control_N.genus, file = "D1_Control_N.genus.csv")

#Day 1 BB nasal total Genus
D1_BB_N.genus <- filter(fobar2.gather, grepl(paste("D1_BB"), All))
D1_BB_N.genus
D1_BB_N.genus$Indicator<-"Zero"
D1_BB_N.genus$Indicator[D1_BB_N.genus$value>0]<-D1_BB_N.genus$Genus[D1_BB_N.genus$value>0]
unique(D1_BB_N.genus$Indicator) #87 for indicator column
write.csv(D1_BB_N.genus, file = "D1_BB_N.genus.csv")

#Day 3 Control nasal total Genus
D3_Control_N.genus <- filter(fobar2.gather, grepl(paste("D3_Control"), All))
D3_Control_N.genus
D3_Control_N.genus$Indicator<-"Zero"
D3_Control_N.genus$Indicator[D3_Control_N.genus$value>0]<-D3_Control_N.genus$Genus[D3_Control_N.genus$value>0]
unique(D3_Control_N.genus$Indicator) #90 for indicator column
write.csv(D3_Control_N.genus, file = "D3_Control_N.genus.csv")

#Day 3 BB nasal total Genus
D3_BB_N.genus <- filter(fobar2.gather, grepl(paste("D3_BB"), All))
D3_BB_N.genus
D3_BB_N.genus$Indicator<-"Zero"
D3_BB_N.genus$Indicator[D3_BB_N.genus$value>0]<-D3_BB_N.genus$Genus[D3_BB_N.genus$value>0]
unique(D3_BB_N.genus$Indicator) #91 for indicator column
write.csv(D3_BB_N.genus, file = "D3_BB_N.genus.csv")

#Day 7 Control nasal total Genus
D7_Control_N.genus <- filter(fobar2.gather, grepl(paste("D7_Control"), All))
D7_Control_N.genus
D7_Control_N.genus$Indicator<-"Zero"
D7_Control_N.genus$Indicator[D7_Control_N.genus$value>0]<-D7_Control_N.genus$Genus[D7_Control_N.genus$value>0]
unique(D7_Control_N.genus$Indicator) #88 for indicator column
write.csv(D7_Control_N.genus, file = "D7_Control_N.genus.csv")

#Day 7 BB nasal total Genus
D7_BB_N.genus <- filter(fobar2.gather, grepl(paste("D7_BB"), All))
D7_BB_N.genus
D7_BB_N.genus$Indicator<-"Zero"
D7_BB_N.genus$Indicator[D7_BB_N.genus$value>0]<-D7_BB_N.genus$Genus[D7_BB_N.genus$value>0]
unique(D7_BB_N.genus$Indicator) #86 for indicator column
write.csv(D7_BB_N.genus, file = "D7_BB_N.genus.csv")

#Day 10 Control nasal total Genus
D10_Control_N.genus <- filter(fobar2.gather, grepl(paste("D10_Control"), All))
D10_Control_N.genus
D10_Control_N.genus$Indicator<-"Zero"
D10_Control_N.genus$Indicator[D10_Control_N.genus$value>0]<-D10_Control_N.genus$Genus[D10_Control_N.genus$value>0]
unique(D10_Control_N.genus$Indicator) #92 for indicator column
write.csv(D10_Control_N.genus, file = "D10_Control_N.genus.csv")

#Day 10 BB nasal total Genus
D10_BB_N.genus <- filter(fobar2.gather, grepl(paste("D10_BB"), All))
D10_BB_N.genus
D10_BB_N.genus$Indicator<-"Zero"
D10_BB_N.genus$Indicator[D10_BB_N.genus$value>0]<-D10_BB_N.genus$Genus[D10_BB_N.genus$value>0]
unique(D10_BB_N.genus$Indicator) #92 for indicator column
write.csv(D10_BB_N.genus, file = "D10_BB_N.genus.csv")

#Day 14 Control nasal total Genus
D14_Control_N.genus <- filter(fobar2.gather, grepl(paste("D14_Control"), All))
D14_Control_N.genus
D14_Control_N.genus$Indicator<-"Zero"
D14_Control_N.genus$Indicator[D14_Control_N.genus$value>0]<-D14_Control_N.genus$Genus[D14_Control_N.genus$value>0]
unique(D14_Control_N.genus$Indicator) #91 for indicator column
write.csv(D14_Control_N.genus, file = "D14_Control_N.genus.csv")

#Day 14 BB nasal total Genus
D14_BB_N.genus <- filter(fobar2.gather, grepl(paste("D14_BB"), All))
D14_BB_N.genus
D14_BB_N.genus$Indicator<-"Zero"
D14_BB_N.genus$Indicator[D14_BB_N.genus$value>0]<-D14_BB_N.genus$Genus[D14_BB_N.genus$value>0]
unique(D14_BB_N.genus$Indicator) #93 for indicator column
write.csv(D14_BB_N.genus, file = "D14_BB_N.genus.csv")

#Day 21 Control nasal total Genus
D21_Control_N.genus <- filter(fobar2.gather, grepl(paste("D21_Control"), All))
D21_Control_N.genus
D21_Control_N.genus$Indicator<-"Zero"
D21_Control_N.genus$Indicator[D21_Control_N.genus$value>0]<-D21_Control_N.genus$Genus[D21_Control_N.genus$value>0]
unique(D21_Control_N.genus$Indicator) #94 for indicator column
write.csv(D21_Control_N.genus, file = "D21_Control_N.genus.csv")

#Day 21 BB nasal total Genus
D21_BB_N.genus <- filter(fobar2.gather, grepl(paste("D21_BB"), All))
D21_BB_N.genus
D21_BB_N.genus$Indicator<-"Zero"
D21_BB_N.genus$Indicator[D21_BB_N.genus$value>0]<-D21_BB_N.genus$Genus[D21_BB_N.genus$value>0]
unique(D21_BB_N.genus$Indicator) #94 for indicator column
write.csv(D21_BB_N.genus, file = "D21_BB_N.genus.csv")

#Day 28 Control nasal total Genus
D28_Control_N.genus <- filter(fobar2.gather, grepl(paste("D28_Control"), All))
D28_Control_N.genus
D28_Control_N.genus$Indicator<-"Zero"
D28_Control_N.genus$Indicator[D28_Control_N.genus$value>0]<-D28_Control_N.genus$Genus[D28_Control_N.genus$value>0]
unique(D28_Control_N.genus$Indicator) #91 for indicator column
write.csv(D28_Control_N.genus, file = "D28_Control_N.genus.csv")

#Day 28 BB nasal total Genus
D28_BB_N.genus <- filter(fobar2.gather, grepl(paste("D28_BB"), All))
D28_BB_N.genus
D28_BB_N.genus$Indicator<-"Zero"
D28_BB_N.genus$Indicator[D28_BB_N.genus$value>0]<-D28_BB_N.genus$Genus[D28_BB_N.genus$value>0]
unique(D28_BB_N.genus$Indicator) #94 for indicator column
write.csv(D28_BB_N.genus, file = "D28_BB_N.genus.csv")

#Day 36 Control nasal total Genus
D36_Control_N.genus <- filter(fobar2.gather, grepl(paste("D36_Control"), All))
D36_Control_N.genus
D36_Control_N.genus$Indicator<-"Zero"
D36_Control_N.genus$Indicator[D36_Control_N.genus$value>0]<-D36_Control_N.genus$Genus[D36_Control_N.genus$value>0]
unique(D36_Control_N.genus$Indicator) #94 for indicator column
write.csv(D36_Control_N.genus, file = "D36_Control_N.genus.csv")

#Day 36 BB nasal total Genus
D36_BB_N.genus <- filter(fobar2.gather, grepl(paste("D36_BB"), All))
D36_BB_N.genus
D36_BB_N.genus$Indicator<-"Zero"
D36_BB_N.genus$Indicator[D36_BB_N.genus$value>0]<-D36_BB_N.genus$Genus[D36_BB_N.genus$value>0]
unique(D36_BB_N.genus$Indicator) #95 for indicator column
write.csv(D36_BB_N.genus, file = "D36_BB_N.genus.csv")

#Day 42 Control nasal total Genus
D42_Control_N.genus <- filter(fobar2.gather, grepl(paste("D42_Control"), All))
D42_Control_N.genus
D42_Control_N.genus$Indicator<-"Zero"
D42_Control_N.genus$Indicator[D42_Control_N.genus$value>0]<-D42_Control_N.genus$Genus[D42_Control_N.genus$value>0]
unique(D42_Control_N.genus$Indicator) #92 for indicator column
write.csv(D42_Control_N.genus, file = "D42_Control_N.genus.csv")

#Day 42 BB nasal total Genus
D42_BB_N.genus <- filter(fobar2.gather, grepl(paste("D42_BB"), All))
D42_BB_N.genus
D42_BB_N.genus$Indicator<-"Zero"
D42_BB_N.genus$Indicator[D42_BB_N.genus$value>0]<-D42_BB_N.genus$Genus[D42_BB_N.genus$value>0]
unique(D42_BB_N.genus$Indicator) #95 for indicator column
write.csv(D42_BB_N.genus, file = "D42_BB_N.genus.csv")


#####Nasal genus DESeq2 abundance######
#Modify nasalgen2.1 to narrow list of genera to the list of DESeq2 significant genera
write.csv(nasalgen2.1, file="SRD129_Nasal_BBControlNoDNEG_GenusAbundanceMatchDESeq2List.csv")
nasalgen2.1.a <- read.csv("SRD129_Nasal_BBControlNoDNEG_GenusAbundanceMatchDESeq2ListModified.csv")
nasalgen2.1.a$Day = factor(nasalgen2.1.a$Day, levels = c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
unique(nasalgen2.1$Treatment) #BB Control

#Nasal genus plot BB and control, abundance plots for specific DESeq2 genera
(nasalgen2.1.aplot <- ggplot(data=nasalgen2.1.a, aes(x=Day, y=Percent.abundance, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, scales = "free", nrow=2) + ylab('Relative abundance (%) at nasal site') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    scale_fill_igv(name = "Genus") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual("Genus", values=c("Control" = "#FF3300", "BB" = "#F8766D")))
nasalgen2.1.aplot