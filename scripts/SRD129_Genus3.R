############################################################
#SRD129 total genera - BB, Control; no DNEG12, DNEG6
#Kathy Mou

#For generating total genera found in each treatment group per day and plot as bar graphs
#Also assessing total number of unique families in each day+treatment group per tissue type
#Only analyzed data without bad samples

#Set working directory
setwd("~/Desktop/Microbiome/Projects/SRD129/SRD129_2000singletons")
setwd("~/Desktop/SRD129/SRD129_2000singletons")
setwd("C:/Users/Kathy.Mou/Desktop/SRD129/SRD129_2000singletons")

#Clear workspace and load necessary packages
rm(list=ls())

#Load libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(cowplot)
library("ggsci")

#Load saved image
load("SRD129.Genus3.RData")



#Save image
save.image(file="SRD129.Genus3.RData")

############################################################

##No bad samples metadata

otu <- import_mothur(mothur_shared_file = 'SRD129.outsingletons.abund.opti_mcc.0.03.subsample.shared')
sample_sums(otu)
taxo <- import_mothur(mothur_constaxonomy_file = 'SRD129.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
shared <- read.table('SRD129.outsingletons.abund.opti_mcc.0.03.subsample.shared', header = TRUE)
meta <- read.table('SRD129metadata06292018nobad.csv', header = TRUE, sep = ",")
head(meta)

#Checking # of samples in each group/day
meta2 <- meta
meta2 <- unite(meta2, "All_Tissue", c("All", "Tissue"), remove=FALSE)
meta2_table <- table(meta2$All_Tissue)
write.csv(meta2_table, "SRD129metadataNoBadSampleCount.csv")
meta3_table <- read.table('SRD129metadataNoBadSampleCount.csv', header = TRUE, sep= ",")
meta3_table <- separate(meta3_table, "Var1", c("Day", "Treatment", "Tissue"), "_")
write.csv(meta3_table, "SRD129metadataNoBadSampleCount.csv")

colnames(meta)[1] <- 'group' #rename as group temporarily. Will use to set as rownames later and remove the
#"group" column
#meta$Day<- gsub("D", "", meta$Day)
meta$group <- as.character(meta$group)
head(meta)
#meta$pignum <- gsub('P([0-9]+)([A-Z]+)D([0-9]+)', '\\1', meta$group)
#meta$tissue <- gsub('P([0-9]+)([A-Z]+)D([0-9]+)', '\\2', meta$group)
#meta$day <- gsub('P([0-9]+)([A-Z]+)D([0-9]+)', '\\3', meta$group)

#meta$pignum <- gsub("SRD129", "", meta$pignum) #remove SRD129 from pignum column
#meta$tissue <- gsub("SRD129", "", meta$tissue) #remove SRD129 from tissue column
#meta$day <- gsub("SRD129", "", meta$day) #remove SRD129 from day column

phy_meta <- sample_data(meta) 
rownames(phy_meta) <- phy_meta$group #Set group names as rownames
head(phy_meta)
phy_meta <- phy_meta[,-1] #remove group column
head(phy_meta)
sample_sums(otu)

SRD129 <- phyloseq(otu, taxo) #build phyloseq-class objects from their components
SRD129 <- merge_phyloseq(SRD129, phy_meta)  # combines the metadata with this phyloseq object
colnames(tax_table(SRD129)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
#SRD129 <- subset_samples(SRD129, experiment == 'SRD129')
#SRD129 <- prune_samples(sample_sums(SRD129) > 700, SRD129)  # This removes samples that have fewer than 700 sequences associated with them.
#700 is arbitrary but you won't need to run this as you already rarified/subsampled your data

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
head(meta)
meta <- meta[,-1]
head(meta)
names(meta)[1] <- paste("group")
head(meta)
fobar <- merge(meta, phyla_tab2, by = 'group') #merge meta with phyla_tab2 by "group"
head(fobar)

fobar.gather <- fobar %>% gather(Genus, value, -(group:All))  #converts to long form dataframe, this is handy for using ggplot faceting functions, check out tidyverse tutorials
#create new columns Genus, value; add columns group to All before Genus and value
head(fobar.gather)

#modified ggplot
ggplot(data=fobar.gather, aes(x=Tissue, y=value, fill=Genus)) +
  geom_bar(stat = 'identity') +
  facet_grid(Treatment~Day) + ylab('Percent of total genus in nasal and tonsil tissue')

#check if there is an extra "group" column. If so, run the next command and modify which column to remove
which(colnames(phyla_tab2)=="group") #result says column 429
phyla_tab3 <- phyla_tab2[,-429] ##drop the 237th column because it's an extra group (sample name) column
phyla_tab4 <- phyla_tab3[,colSums(phyla_tab3)>0.1] #keep the columns that have greater than 0.1 value
phyla_tab4$group <- rownames(phyla_tab4) #rename rownames as group
fobar2 <- merge(meta, phyla_tab4, by = 'group')
head(fobar2)

#To see how many T are in meta$Tissue: 
length(which(meta$Tissue== "T")) #431

#To see how many N are in meta$Tissue:
length(which(meta$Tissue== "N")) #460

fobar2.gather <- fobar2 %>% gather(Genus, value, -(group:All))  #converts to long form dataframe, this is handy for using ggplot faceting functions, check out tidyverse tutorials
#fobar2.gather$day <- as.numeric(fobar2.gather$Day)
head(fobar2.gather)
#fobar2.gather$day <- NULL

#reordering days 0-14 in plot
levels(sample_data(fobar2.gather)$Day) #Output: "D0"     "D1"     "D10"    "D14"    "D21"    "D28"    "D3"     "D36"    "D42"    "D7"     "DNEG12" "DNEG6" 
fobar2.gather$Day <- factor(fobar2.gather$Day, levels=c("DNEG12", "DNEG6", "D0",  "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
head(fobar2.gather$Day) #DNEG12 DNEG6 D0 D1 D3 D7 D10 D14 D21 D28 D36 D42

fobar2.gather$All <- paste(fobar2.gather$Day, fobar2.gather$Treatment, fobar2.gather$Tissue, sep = '_')
#Create "All" column with Day, Treatment and Tissue
#Use this to make % total community plot
head(fobar2.gather$All)

#Ways to summarize unique genera
unique(fobar2.gather$Genus) #Total number of unique genera: 181
#OR
fobar2.gather %>% summarise_each(funs(n_distinct)) #181 total unique genus

fobar2.gather <- fobar2.gather %>% group_by(All) %>% mutate(value2=(value/(length(All)/181))*100)
#181 refers to number of Genus (counted based on summarise_each function that was ran on fobar2.gather
#to get total number unique items under "Genus" column)

theme_update(plot.title = element_text(hjust = 0.5)) #center title
ggplot(data=fobar2.gather, aes(x=Tissue, y=value2, fill=Genus)) +
  geom_bar(stat = 'identity') +
  facet_grid(Treatment~Day) + ylab('Percent of total community') +
  ggtitle("Total percent genus members in nasal wash or tonsil tissue \namong 3 treatment groups")


#############NASAL#############
#Subset each "All" group from fobar2.gather and save as individual csv files.
#For example:
D0_Control_NW.genus <- subset(fobar2.gather, All=="D0_Control_NW") #EDIT THIS
write.csv(D0_Control_NW.genus, file="D0_Control_NW.genus.csv")  #EDIT THIS

#Calculate the total % percent abundance of each genera on each sample (I used JMP to do this) 
#and save results in a spreadsheet editor such as Excel (see D0_Control_NW.genus.xlsx for an example)
#Since we are only interested in genera that are above 2% abundance, 
#calculate total percentage of all other genera that are less than 2% in the spreadsheet and label as "Other".
#Create a new spreadsheet and label as "Nasal genus.csv" or "Tonsil genus.csv". 
#Create the following columns: Day, Treatment group, Tissue, Percent Abundance, and Genus. 
#Copy the list of genera and their percent abundances from each of the individual Excel files to the respective "Nasal genus.csv" or "Tonsil genus.csv" spreadsheet.
#Fill in the other columns manually (Day, Treatment Group, Tissue). 
#You should have a file similar to SRD129_Nasal_GenusPercentAbundance.csv. Continue to the next step.

#Jan. 30, 2019: edit this csv file to include updated D1 BB percentage info
nasalgen <- read.table('SRD129_Nasal_GenusPercentAbundance.csv', header = TRUE, sep = ",")
head(nasalgen)



#####Nasal genus with other category above 2% abundance######
nasalgen$More.than.2=as.character(nasalgen$Genus)
str(nasalgen$More.than.2) #Compactly display the internal structure of an R object,
nasalgen$More.than.2[nasalgen$Percent.abundance<2]<-"Other"
write.csv(nasalgen, file = "SRD129_Nasal_GenusPercentAbundanceAbove2percent.csv")
#Edited this file and saved as "SRD129_Nasal_GenusPercentAbundanceAbove2percentAddedTo100.xlsx"
#Created a new file called SRD129_Nasal_GenusPercentAbundance2.csv, modified column names, and imported as nasalgen2
nasalgen2 <- read.table('SRD129_Nasal_GenusPercentAbundance2.csv', header = TRUE, sep = ",")
head(nasalgen2)
unique(nasalgen2$Day)
unique(nasalgen2$Treatment)
nasalgen2$Day = factor(nasalgen2$Day, levels = c("DNEG12", "DNEG6", "D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
nasalgen2.2<- nasalgen2[!grepl("PRRSV", nasalgen2$Treatment),]
nasalgen2.2<- nasalgen2.2[!grepl("IAV", nasalgen2.2$Treatment),]
nasalgen2.2<- nasalgen2.2[!grepl("DNEG12", nasalgen2.2$Day),]
nasalgen2.2<- nasalgen2.2[!grepl("DNEG6", nasalgen2.2$Day),]
unique(nasalgen2.2$Day) #D0  D1  D3  D7  D10 D14 D21 D28 D36 D42
unique(nasalgen2.2$Treatment) #Control BB

#Nasal genus plot BB and control only, abundance plots for each genus, more than 2%
(nasalgen2.2plot <- ggplot(data=nasalgen2.2, aes(x=Day, y=Percent.Abundance, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, scales = "free") + ylab('Relative abundance (%) at nasal site') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 2)) +
    theme(axis.line = element_line()) +
    scale_fill_igv(name = "Genus") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
nasalgen2.2plot

#Nasal genus plot BB and control only, abundance plot for each day, more than 2%
(nasalgen2.2bplot <- ggplot(data=nasalgen2.2, aes(x=Treatment, y=Percent.Abundance, fill=Genus)) +
    geom_bar(stat = 'identity') +
    #geom_bar(stat= 'identity', colour='black') +
    #theme(legend.key = element_rect = element_rect(colour='black', size=1.5)) +
    facet_grid(~Day) + ylab('Relative abundance (%) at nasal site') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank()) +
    scale_fill_igv(name = "Genus") +
    theme(legend.direction = "vertical"))




#####Nasal genus all genera######
#Nasalgen2.1 (BB and control only, all genera)
nasalgen2.1 <- nasalgen[!grepl("PRRSV", nasalgen$Treatment),]
nasalgen2.1 <- nasalgen2.1[!grepl("IAV", nasalgen2.1$Treatment),]
nasalgen2.1<- nasalgen2.1[!grepl("DNEG12", nasalgen2.1$Day),]
nasalgen2.1<- nasalgen2.1[!grepl("DNEG6", nasalgen2.1$Day),]
unique(nasalgen2.1$Treatment) #Control BB
unique(nasalgen2.1$Day) #D0  D1  D10 D14 D21 D28 D3  D36 D42 D7
nasalgen2.1$Day = factor(nasalgen2.1$Day, levels = c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
write.csv(nasalgen2.1, file = "SRD129_Nasal_BBControlNoDNEGAllGenusPercentAbundance.csv")


#Nasal genus plot BB and control only, abundance plots for each genus, all genera
(nasalgen2.1plot <- ggplot(data=nasalgen2.1, aes(x=Day, y=Percent.abundance, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, scales = "free") + ylab('Relative abundance (%) at nasal site') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 2)) +
    theme(axis.line = element_line()) +
    scale_fill_igv(name = "Genus") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
nasalgen2.1plot

#Nasal genus plot BB and control only, abundance plot for each day, all genera
(nasalgen2.1bplot <- ggplot(data=nasalgen2.1, aes(x=Treatment, y=Percent.abundance, fill=Genus)) +
    geom_bar(stat = 'identity') +
    #geom_bar(stat= 'identity', colour='black') +
    #theme(legend.key = element_rect = element_rect(colour='black', size=1.5)) +
    facet_grid(~Day) + ylab('Relative abundance (%) at nasal site') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank()) +
    theme(legend.direction = "vertical"))
nasalgen2.1bplot





#####Nasal genus with other category above 1% abundance######
#Nasalgen.2 (BB and control only, more than 1% genera)
nasalgen.2<- nasalgen[!grepl("PRRSV", nasalgen$Treatment),]
nasalgen.2<- nasalgen.2[!grepl("IAV", nasalgen.2$Treatment),]
nasalgen.2 <- nasalgen.2[!grepl("DNEG12", nasalgen.2$Day),]
nasalgen.2 <- nasalgen.2[!grepl("DNEG6", nasalgen.2$Day),]
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


#Nasalgen.2 abundance plot for each day, BB and control only, no DNEG12/6, more than 1% genera
(nasalgen.2bplot <- ggplot(data=nasalgen.2b, aes(x=Treatment, y=Percent.abundance, fill=More.than.2)) +
    geom_bar(stat = 'identity') +
    #geom_bar(stat= 'identity', colour='black') +
    #theme(legend.key = element_rect = element_rect(colour='black', size=1.5)) +
    facet_grid(~Day) + ylab('Relative abundance (%) at nasal site') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank()) +
    scale_fill_igv(name = "Genus") +
    theme(legend.direction = "vertical"))
nasalgen.2bplot

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






#####Nasal genus DESeq2 abundance######
#Modify nasalgen2.1 to narrow list of genera to the list of DESeq2 significant genera
write.csv(nasalgen2.1, file="SRD129_Nasal_BBControlNoDNEG_GenusAbundanceMatchDESeq2List.csv")
nasalgen2.1.a <- read.csv("SRD129_Nasal_BBControlNoDNEG_GenusAbundanceMatchDESeq2ListModified.csv")
nasalgen2.1.a$Day = factor(nasalgen2.1.a$Day, levels = c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
unique(nasalgen2.1$Treatment) #BB Control

#Nasal genus plot flu and control only, abundance plots for specific DESeq2 genera
(nasalgen2.1.aplot <- ggplot(data=nasalgen2.1.a, aes(x=Day, y=Percent.abundance, fill=Treatment)) +
    geom_bar(stat = 'identity', position="dodge") +
    facet_wrap(~Genus, scales = "free", nrow=2) + ylab('Relative abundance (%) at nasal site') +
    theme_gray()+
    theme(plot.title = element_text(hjust = 3)) +
    theme(axis.line = element_line()) +
    scale_fill_igv(name = "Genus") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual("Genus", values=c("Control" = "#FF3300", "IAV" = "#00CCCC")))
nasalgen2.1.aplot



############################################## Nasal total Genus ##############################################################

#No bad samples

#get names of all groups in "All" column of fobar2.gather
fobar2.gather$All

#Day -12 Control nasal total genus
DNEG12_Control_N.genus <- filter(fobar2.gather, grepl(paste("DNEG12_Control_N"), All))
DNEG12_Control_N.genus
DNEG12_Control_N.genus$Indicator<-"Zero"
DNEG12_Control_N.genus$Indicator[DNEG12_Control_N.genus$value>0]<-DNEG12_Control_N.genus$Genus[DNEG12_Control_N.genus$value>0]
unique(DNEG12_Control_N.genus$Indicator) #132 for indicator column
DNEG12_Control_N.genus %>% summarise_each(funs(n_distinct)) #132 for Indicator column
write.csv(DNEG12_Control_N.genus, file = "DNEG12_Control_N.genus.csv")

#Day -12 BB nasal total Genus
DNEG12_BB_N.genus <- filter(fobar2.gather, grepl(paste("DNEG12_BB_N"), All))
DNEG12_BB_N.genus
DNEG12_BB_N.genus$Indicator<-"Zero"
DNEG12_BB_N.genus$Indicator[DNEG12_BB_N.genus$value>0]<-DNEG12_BB_N.genus$Genus[DNEG12_BB_N.genus$value>0]
unique(DNEG12_BB_N.genus$Indicator) #148 for indicator column
write.csv(DNEG12_BB_N.genus, file = "DNEG12_BB_N.genus.csv")

#Day -12 IAV nasal total Genus
DNEG12_IAV_N.genus <- filter(fobar2.gather, grepl(paste("DNEG12_IAV_N"), All))
DNEG12_IAV_N.genus
DNEG12_IAV_N.genus$Indicator<-"Zero"
DNEG12_IAV_N.genus$Indicator[DNEG12_IAV_N.genus$value>0]<-DNEG12_IAV_N.genus$Genus[DNEG12_IAV_N.genus$value>0]
unique(DNEG12_IAV_N.genus$Indicator) #138 for indicator column
write.csv(DNEG12_IAV_N.genus, file = "DNEG12_IAV_N.genus.csv")

#Day -12 PRRSV nasal total Genus
DNEG12_PRRSV_N.genus <- filter(fobar2.gather, grepl(paste("DNEG12_PRRSV_N"), All))
DNEG12_PRRSV_N.genus
DNEG12_PRRSV_N.genus$Indicator<-"Zero"
DNEG12_PRRSV_N.genus$Indicator[DNEG12_PRRSV_N.genus$value>0]<-DNEG12_PRRSV_N.genus$Genus[DNEG12_PRRSV_N.genus$value>0]
unique(DNEG12_PRRSV_N.genus$Indicator) #138 for indicator column
write.csv(DNEG12_PRRSV_N.genus, file = "DNEG12_PRRSV_N.genus.csv")

DNEG12.genus <- filter(fobar2.gather, grepl(paste("DNEG12"), Day))
summary(DNEG12.genus$Treatment) #Tonsil and Nasal
#  BB Control     IAV   PRRSV 
#3620    3620    3620    3620

#Making sure nothing got left out when subsetting data
grep("Control",DNEG12_Control_N.genus$Treatment) #1810
grep("Control",DNEG12_Control_T.genus$Treatment) #1810
grep("BB",DNEG12_BB_N.genus$Treatment) #1810
grep("BB",DNEG12_BB_T.genus$Treatment) #1810
grep("IAV",DNEG12_IAV_N.genus$Treatment) #1810
grep("IAV",DNEG12_IAV_T.genus$Treatment) #1810
grep("PRRSV",DNEG12_PRRSV_N.genus$Treatment) #1810
grep("PRRSV",DNEG12_PRRSV_T.genus$Treatment) #1810





#Day -6 Control nasal total genus
DNEG6_Control_N.genus <- filter(fobar2.gather, grepl(paste("DNEG6_Control_N"), All))
DNEG6_Control_N.genus
DNEG6_Control_N.genus$Indicator<-"Zero"
DNEG6_Control_N.genus$Indicator[DNEG6_Control_N.genus$value>0]<-DNEG6_Control_N.genus$Genus[DNEG6_Control_N.genus$value>0]
unique(DNEG6_Control_N.genus$Indicator) #148 for indicator column
write.csv(DNEG6_Control_N.genus, file = "DNEG6_Control_N.genus.csv")

#Day -6 BB nasal total Genus
DNEG6_BB_N.genus <- filter(fobar2.gather, grepl(paste("DNEG6_BB_N"), All))
DNEG6_BB_N.genus
DNEG6_BB_N.genus$Indicator<-"Zero"
DNEG6_BB_N.genus$Indicator[DNEG6_BB_N.genus$value>0]<-DNEG6_BB_N.genus$Genus[DNEG6_BB_N.genus$value>0]
unique(DNEG6_BB_N.genus$Indicator) #153 for indicator column
write.csv(DNEG6_BB_N.genus, file = "DNEG6_BB_N.genus.csv")

#Day -6 IAV nasal total Genus
DNEG6_IAV_N.genus <- filter(fobar2.gather, grepl(paste("DNEG6_IAV_N"), All))
DNEG6_IAV_N.genus
DNEG6_IAV_N.genus$Indicator<-"Zero"
DNEG6_IAV_N.genus$Indicator[DNEG6_IAV_N.genus$value>0]<-DNEG6_IAV_N.genus$Genus[DNEG6_IAV_N.genus$value>0]
unique(DNEG6_IAV_N.genus$Indicator) #145 for indicator column
write.csv(DNEG6_IAV_N.genus, file = "DNEG6_IAV_N.genus.csv")

#Day -6 PRRSV nasal total Genus
DNEG6_PRRSV_N.genus <- filter(fobar2.gather, grepl(paste("DNEG6_PRRSV_N"), All))
DNEG6_PRRSV_N.genus
DNEG6_PRRSV_N.genus$Indicator<-"Zero"
DNEG6_PRRSV_N.genus$Indicator[DNEG6_PRRSV_N.genus$value>0]<-DNEG6_PRRSV_N.genus$Genus[DNEG6_PRRSV_N.genus$value>0]
unique(DNEG6_PRRSV_N.genus$Indicator) #153 for indicator column
write.csv(DNEG6_PRRSV_N.genus, file = "DNEG6_PRRSV_N.genus.csv")

DNEG6.genus <- filter(fobar2.gather, grepl(paste("DNEG6"), Day))
summary(DNEG6.genus$Treatment) #Tonsil and Nasal
#  BB Control     IAV   PRRSV 
#2534    2715    2715    2172

#Making sure nothing got left out when subsetting data
grep("Control",DNEG6_Control_N.genus$Treatment) #1810
grep("Control",DNEG6_Control_T.genus$Treatment) #9.5
grep("BB",DNEG6_BB_N.genus$Treatment) #1629
grep("BB",DNEG6_BB_T.genus$Treatment) #905
grep("IAV",DNEG6_IAV_N.genus$Treatment) #1810
grep("IAV",DNEG6_IAV_T.genus$Treatment) #905
grep("PRRSV",DNEG6_PRRSV_N.genus$Treatment) #1810
grep("PRRSV",DNEG6_PRRSV_T.genus$Treatment) #362





fobar2.gather$All
#Day 0 Control nasal total Genus
D0_Control_N.genus <- filter(fobar2.gather, grepl(paste("D0_Control_N"), All))
D0_Control_N.genus
D0_Control_N.genus$Indicator<-"Zero"
D0_Control_N.genus$Indicator[D0_Control_N.genus$value>0]<-D0_Control_N.genus$Genus[D0_Control_N.genus$value>0]
unique(D0_Control_N.genus$Indicator) #148 for indicator column
write.csv(D0_Control_N.genus, file = "D0_Control_N.genus.csv")

#Day 0 BB nasal total Genus
D0_BB_N.genus <- filter(fobar2.gather, grepl(paste("D0_BB_N"), All))
D0_BB_N.genus
D0_BB_N.genus$Indicator<-"Zero"
D0_BB_N.genus$Indicator[D0_BB_N.genus$value>0]<-D0_BB_N.genus$Genus[D0_BB_N.genus$value>0]
unique(D0_BB_N.genus$Indicator) #134 for indicator column
write.csv(D0_BB_N.genus, file = "D0_BB_N.genus.csv")

#Day 0 IAV nasal total Genus
D0_IAV_N.genus <- filter(fobar2.gather, grepl(paste("D0_IAV_N"), All))
D0_IAV_N.genus
D0_IAV_N.genus$Indicator<-"Zero"
D0_IAV_N.genus$Indicator[D0_IAV_N.genus$value>0]<-D0_IAV_N.genus$Genus[D0_IAV_N.genus$value>0]
unique(D0_IAV_N.genus$Indicator) #135 for indicator column
write.csv(D0_IAV_N.genus, file = "D0_IAV_N.genus.csv")

#Day 0 PRRSV nasal total Genus
D0_PRRSV_N.genus <- filter(fobar2.gather, grepl(paste("D0_PRRSV_N"), All))
D0_PRRSV_N.genus
D0_PRRSV_N.genus$Indicator<-"Zero"
D0_PRRSV_N.genus$Indicator[D0_PRRSV_N.genus$value>0]<-D0_PRRSV_N.genus$Genus[D0_PRRSV_N.genus$value>0]
unique(D0_PRRSV_N.genus$Indicator) #129 for indicator column
write.csv(D0_PRRSV_N.genus, file = "D0_PRRSV_N.genus.csv")

D0.genus <- filter(fobar2.gather, grepl(paste("D0"), Day))
summary(D0.genus$Treatment) #Tonsil and Nasal
#  BB Control     IAV   PRRSV 
#3439    3439    3258    3258 






#Day 1 Control nasal total Genus
D1_Control_N.genus <- filter(fobar2.gather, grepl(paste("D1_Control_N"), All))
D1_Control_N.genus
D1_Control_N.genus$Indicator<-"Zero"
D1_Control_N.genus$Indicator[D1_Control_N.genus$value>0]<-D1_Control_N.genus$Genus[D1_Control_N.genus$value>0]
unique(D1_Control_N.genus$Indicator) #140 for indicator column
write.csv(D1_Control_N.genus, file = "D1_Control_N.genus.csv")

#Day 1 BB nasal total Genus
D1_BB_N.genus <- filter(fobar2.gather, grepl(paste("D1_BB_N"), All))
D1_BB_N.genus
D1_BB_N.genus$Indicator<-"Zero"
D1_BB_N.genus$Indicator[D1_BB_N.genus$value>0]<-D1_BB_N.genus$Genus[D1_BB_N.genus$value>0]
unique(D1_BB_N.genus$Indicator) #136 for indicator column
write.csv(D1_BB_N.genus, file = "D1_BB_N.genus.csv")
#Fixing the 200%
D1_BB_N2.genus<- D1_BB_N.genus[!grepl("NT", D1_BB_N.genus$Tissue),]
D1_BB_N2.genus$Indicator[D1_BB_N2.genus$value>0]<-D1_BB_N2.genus$Genus[D1_BB_N2.genus$value>0]
unique(D1_BB_N2.genus$Indicator) #133 for indicator column
write.csv(D1_BB_N2.genus, file = "D1_BB_N2.genus.csv")

#Day 1 IAV nasal total Genus
D1_IAV_N.genus <- filter(fobar2.gather, grepl(paste("D1_IAV_N"), All))
D1_IAV_N.genus
D1_IAV_N.genus$Indicator<-"Zero"
D1_IAV_N.genus$Indicator[D1_IAV_N.genus$value>0]<-D1_IAV_N.genus$Genus[D1_IAV_N.genus$value>0]
unique(D1_IAV_N.genus$Indicator) #137 for indicator column
write.csv(D1_IAV_N.genus, file = "D1_IAV_N.genus.csv")

#Day 1 PRRSV nasal total Genus
D1_PRRSV_N.genus <- filter(fobar2.gather, grepl(paste("D1_PRRSV_N"), All))
D1_PRRSV_N.genus
D1_PRRSV_N.genus$Indicator<-"Zero"
D1_PRRSV_N.genus$Indicator[D1_PRRSV_N.genus$value>0]<-D1_PRRSV_N.genus$Genus[D1_PRRSV_N.genus$value>0]
unique(D1_PRRSV_N.genus$Indicator) #136 for indicator column
write.csv(D1_PRRSV_N.genus, file = "D1_PRRSV_N.genus.csv")

D1.genus <- filter(fobar2.gather, grepl(paste("D1_(*)"), All))
summary(D1.genus$Treatment) #Tonsil and Nasal
#   BB Control     IAV   PRRSV 
# 3439    3620    3439    3439 




#Day 3 Control nasal total Genus
D3_Control_N.genus <- filter(fobar2.gather, grepl(paste("D3_Control_N"), All))
D3_Control_N.genus
D3_Control_N.genus$Indicator<-"Zero"
D3_Control_N.genus$Indicator[D3_Control_N.genus$value>0]<-D3_Control_N.genus$Genus[D3_Control_N.genus$value>0]
unique(D3_Control_N.genus$Indicator) #150 for indicator column
write.csv(D3_Control_N.genus, file = "D3_Control_N.genus.csv")

#Day 3 BB nasal total Genus
D3_BB_N.genus <- filter(fobar2.gather, grepl(paste("D3_BB_N"), All))
D3_BB_N.genus
D3_BB_N.genus$Indicator<-"Zero"
D3_BB_N.genus$Indicator[D3_BB_N.genus$value>0]<-D3_BB_N.genus$Genus[D3_BB_N.genus$value>0]
unique(D3_BB_N.genus$Indicator) #140 for indicator column
write.csv(D3_BB_N.genus, file = "D3_BB_N.genus.csv")

#Day 3 IAV nasal total Genus
D3_IAV_N.genus <- filter(fobar2.gather, grepl(paste("D3_IAV_N"), All))
D3_IAV_N.genus
D3_IAV_N.genus$Indicator<-"Zero"
D3_IAV_N.genus$Indicator[D3_IAV_N.genus$value>0]<-D3_IAV_N.genus$Genus[D3_IAV_N.genus$value>0]
unique(D3_IAV_N.genus$Indicator) #139 for indicator column
write.csv(D3_IAV_N.genus, file = "D3_IAV_N.genus.csv")

#Day 3 PRRSV nasal total Genus
D3_PRRSV_N.genus <- filter(fobar2.gather, grepl(paste("D3_PRRSV_N"), All))
D3_PRRSV_N.genus
D3_PRRSV_N.genus$Indicator<-"Zero"
D3_PRRSV_N.genus$Indicator[D3_PRRSV_N.genus$value>0]<-D3_PRRSV_N.genus$Genus[D3_PRRSV_N.genus$value>0]
unique(D3_PRRSV_N.genus$Indicator) #133 for indicator column
write.csv(D3_PRRSV_N.genus, file = "D3_PRRSV_N.genus.csv")

D3.genus <- filter(fobar2.gather, grepl(paste("D3_(*)"), All))
summary(D3.genus$Treatment) #Tonsil and Nasal
#  BB Control     IAV   PRRSV 
#3439    3620    3620    3620 





#Day 7 Control nasal total Genus
D7_Control_N.genus <- filter(fobar2.gather, grepl(paste("D7_Control_N"), All))
D7_Control_N.genus
D7_Control_N.genus$Indicator<-"Zero"
D7_Control_N.genus$Indicator[D7_Control_N.genus$value>0]<-D7_Control_N.genus$Genus[D7_Control_N.genus$value>0]
unique(D7_Control_N.genus$Indicator) #135 for indicator column
write.csv(D7_Control_N.genus, file = "D7_Control_N.genus.csv")

#Day 7 BB nasal total Genus
D7_BB_N.genus <- filter(fobar2.gather, grepl(paste("D7_BB_N"), All))
D7_BB_N.genus
D7_BB_N.genus$Indicator<-"Zero"
D7_BB_N.genus$Indicator[D7_BB_N.genus$value>0]<-D7_BB_N.genus$Genus[D7_BB_N.genus$value>0]
unique(D7_BB_N.genus$Indicator) #145 for indicator column
write.csv(D7_BB_N.genus, file = "D7_BB_N.genus.csv")

#Day 7 IAV nasal total Genus
D7_IAV_N.genus <- filter(fobar2.gather, grepl(paste("D7_IAV_N"), All))
D7_IAV_N.genus
D7_IAV_N.genus$Indicator<-"Zero"
D7_IAV_N.genus$Indicator[D7_IAV_N.genus$value>0]<-D7_IAV_N.genus$Genus[D7_IAV_N.genus$value>0]
unique(D7_IAV_N.genus$Indicator) #115 for indicator column
write.csv(D7_IAV_N.genus, file = "D7_IAV_N.genus.csv")

#Day 7 PRRSV nasal total Genus
D7_PRRSV_N.genus <- filter(fobar2.gather, grepl(paste("D7_PRRSV_N"), All))
D7_PRRSV_N.genus
D7_PRRSV_N.genus$Indicator<-"Zero"
D7_PRRSV_N.genus$Indicator[D7_PRRSV_N.genus$value>0]<-D7_PRRSV_N.genus$Genus[D7_PRRSV_N.genus$value>0]
unique(D7_PRRSV_N.genus$Indicator) #122 for indicator column
write.csv(D7_PRRSV_N.genus, file = "D7_PRRSV_N.genus.csv")

D7.genus <- filter(fobar2.gather, grepl(paste("D7"), Day))
summary(D7.genus$Treatment) #Tonsil and Nasal
#  BB Control     IAV   PRRSV 
#3439    3620    3620    3620




#Day 10 Control nasal total Genus
D10_Control_N.genus <- filter(fobar2.gather, grepl(paste("D10_Control_N"), All))
D10_Control_N.genus
D10_Control_N.genus$Indicator<-"Zero"
D10_Control_N.genus$Indicator[D10_Control_N.genus$value>0]<-D10_Control_N.genus$Genus[D10_Control_N.genus$value>0]
unique(D10_Control_N.genus$Indicator) #143 for indicator column
write.csv(D10_Control_N.genus, file = "D10_Control_N.genus.csv")

#Day 10 BB nasal total Genus
D10_BB_N.genus <- filter(fobar2.gather, grepl(paste("D10_BB_N"), All))
D10_BB_N.genus
D10_BB_N.genus$Indicator<-"Zero"
D10_BB_N.genus$Indicator[D10_BB_N.genus$value>0]<-D10_BB_N.genus$Genus[D10_BB_N.genus$value>0]
unique(D10_BB_N.genus$Indicator) #143 for indicator column
write.csv(D10_BB_N.genus, file = "D10_BB_N.genus.csv")

#Day 10 IAV nasal total Genus
D10_IAV_N.genus <- filter(fobar2.gather, grepl(paste("D10_IAV_N"), All))
D10_IAV_N.genus
D10_IAV_N.genus$Indicator<-"Zero"
D10_IAV_N.genus$Indicator[D10_IAV_N.genus$value>0]<-D10_IAV_N.genus$Genus[D10_IAV_N.genus$value>0]
unique(D10_IAV_N.genus$Indicator) #139 for indicator column
write.csv(D10_IAV_N.genus, file = "D10_IAV_N.genus.csv")

#Day 10 PRRSV nasal total Genus
D10_PRRSV_N.genus <- filter(fobar2.gather, grepl(paste("D10_PRRSV_N"), All))
D10_PRRSV_N.genus
D10_PRRSV_N.genus$Indicator<-"Zero"
D10_PRRSV_N.genus$Indicator[D10_PRRSV_N.genus$value>0]<-D10_PRRSV_N.genus$Genus[D10_PRRSV_N.genus$value>0]
unique(D10_PRRSV_N.genus$Indicator) #148 for indicator column
write.csv(D10_PRRSV_N.genus, file = "D10_PRRSV_N.genus.csv")

D10.genus <- filter(fobar2.gather, grepl(paste("D10"), Day))
summary(D10.genus$Treatment) #Tonsil and Nasal
#  BB Control     IAV   PRRSV 
#3620    3620    2896    3439




fobar2.gather$All
#Day 14 Control nasal total Genus
D14_Control_N.genus <- filter(fobar2.gather, grepl(paste("D14_Control_N"), All))
D14_Control_N.genus
D14_Control_N.genus$Indicator<-"Zero"
D14_Control_N.genus$Indicator[D14_Control_N.genus$value>0]<-D14_Control_N.genus$Genus[D14_Control_N.genus$value>0]
unique(D14_Control_N.genus$Indicator) #151 for indicator column
write.csv(D14_Control_N.genus, file = "D14_Control_N.genus.csv")

#Day 14 BB nasal total Genus
D14_BB_N.genus <- filter(fobar2.gather, grepl(paste("D14_BB_N"), All))
D14_BB_N.genus
D14_BB_N.genus$Indicator<-"Zero"
D14_BB_N.genus$Indicator[D14_BB_N.genus$value>0]<-D14_BB_N.genus$Genus[D14_BB_N.genus$value>0]
unique(D14_BB_N.genus$Indicator) #149 for indicator column
write.csv(D14_BB_N.genus, file = "D14_BB_N.genus.csv")

#Day 14 IAV nasal total Genus
D14_IAV_N.genus <- filter(fobar2.gather, grepl(paste("D14_IAV_N"), All))
D14_IAV_N.genus
D14_IAV_N.genus$Indicator<-"Zero"
D14_IAV_N.genus$Indicator[D14_IAV_N.genus$value>0]<-D14_IAV_N.genus$Genus[D14_IAV_N.genus$value>0]
unique(D14_IAV_N.genus$Indicator) #155 for indicator column
write.csv(D14_IAV_N.genus, file = "D14_IAV_N.genus.csv")

#Day 14 PRRSV nasal total Genus
D14_PRRSV_N.genus <- filter(fobar2.gather, grepl(paste("D14_PRRSV_N"), All))
D14_PRRSV_N.genus
D14_PRRSV_N.genus$Indicator<-"Zero"
D14_PRRSV_N.genus$Indicator[D14_PRRSV_N.genus$value>0]<-D14_PRRSV_N.genus$Genus[D14_PRRSV_N.genus$value>0]
unique(D14_PRRSV_N.genus$Indicator) #144 for indicator column
write.csv(D14_PRRSV_N.genus, file = "D14_PRRSV_N.genus.csv")

D14.genus <- filter(fobar2.gather, grepl(paste("D14"), Day))
summary(D14.genus$Treatment) #Tonsil and Nasal
#  BB Control     IAV   PRRSV 
#3620    3620    3620    3620






#Day 21 Control nasal total Genus
D21_Control_N.genus <- filter(fobar2.gather, grepl(paste("D21_Control_N"), All))
D21_Control_N.genus
D21_Control_N.genus$Indicator<-"Zero"
D21_Control_N.genus$Indicator[D21_Control_N.genus$value>0]<-D21_Control_N.genus$Genus[D21_Control_N.genus$value>0]
unique(D21_Control_N.genus$Indicator) #151 for indicator column
write.csv(D21_Control_N.genus, file = "D21_Control_N.genus.csv")

#Day 21 BB nasal total Genus
D21_BB_N.genus <- filter(fobar2.gather, grepl(paste("D21_BB_N"), All))
D21_BB_N.genus
D21_BB_N.genus$Indicator<-"Zero"
D21_BB_N.genus$Indicator[D21_BB_N.genus$value>0]<-D21_BB_N.genus$Genus[D21_BB_N.genus$value>0]
unique(D21_BB_N.genus$Indicator) #155 for indicator column
write.csv(D21_BB_N.genus, file = "D21_BB_N.genus.csv")

#Day 21 IAV nasal total Genus
D21_IAV_N.genus <- filter(fobar2.gather, grepl(paste("D21_IAV_N"), All))
D21_IAV_N.genus
D21_IAV_N.genus$Indicator<-"Zero"
D21_IAV_N.genus$Indicator[D21_IAV_N.genus$value>0]<-D21_IAV_N.genus$Genus[D21_IAV_N.genus$value>0]
unique(D21_IAV_N.genus$Indicator) #156 for indicator column
write.csv(D21_IAV_N.genus, file = "D21_IAV_N.genus.csv")

#Day 21 PRRSV nasal total Genus
D21_PRRSV_N.genus <- filter(fobar2.gather, grepl(paste("D21_PRRSV_N"), All))
D21_PRRSV_N.genus
D21_PRRSV_N.genus$Indicator<-"Zero"
D21_PRRSV_N.genus$Indicator[D21_PRRSV_N.genus$value>0]<-D21_PRRSV_N.genus$Genus[D21_PRRSV_N.genus$value>0]
unique(D21_PRRSV_N.genus$Indicator) #150 for indicator column
write.csv(D21_PRRSV_N.genus, file = "D21_PRRSV_N.genus.csv")

D21.genus <- filter(fobar2.gather, grepl(paste("D21"), Day))
summary(D21.genus$Treatment) #Tonsil and Nasal
#   BB Control     IAV   PRRSV 
# 3620    3439    2896    3258 






#Day 28 Control nasal total Genus
D28_Control_N.genus <- filter(fobar2.gather, grepl(paste("D28_Control_N"), All))
D28_Control_N.genus
D28_Control_N.genus$Indicator<-"Zero"
D28_Control_N.genus$Indicator[D28_Control_N.genus$value>0]<-D28_Control_N.genus$Genus[D28_Control_N.genus$value>0]
unique(D28_Control_N.genus$Indicator) #139 for indicator column
write.csv(D28_Control_N.genus, file = "D28_Control_N.genus.csv")

#Day 28 BB nasal total Genus
D28_BB_N.genus <- filter(fobar2.gather, grepl(paste("D28_BB_N"), All))
D28_BB_N.genus
D28_BB_N.genus$Indicator<-"Zero"
D28_BB_N.genus$Indicator[D28_BB_N.genus$value>0]<-D28_BB_N.genus$Genus[D28_BB_N.genus$value>0]
unique(D28_BB_N.genus$Indicator) #147 for indicator column
write.csv(D28_BB_N.genus, file = "D28_BB_N.genus.csv")

#Day 28 IAV nasal total Genus
D28_IAV_N.genus <- filter(fobar2.gather, grepl(paste("D28_IAV_N"), All))
D28_IAV_N.genus
D28_IAV_N.genus$Indicator<-"Zero"
D28_IAV_N.genus$Indicator[D28_IAV_N.genus$value>0]<-D28_IAV_N.genus$Genus[D28_IAV_N.genus$value>0]
unique(D28_IAV_N.genus$Indicator) #137 for indicator column
write.csv(D28_IAV_N.genus, file = "D28_IAV_N.genus.csv")

#Day 28 PRRSV nasal total Genus
D28_PRRSV_N.genus <- filter(fobar2.gather, grepl(paste("D28_PRRSV_N"), All))
D28_PRRSV_N.genus
D28_PRRSV_N.genus$Indicator<-"Zero"
D28_PRRSV_N.genus$Indicator[D28_PRRSV_N.genus$value>0]<-D28_PRRSV_N.genus$Genus[D28_PRRSV_N.genus$value>0]
unique(D28_PRRSV_N.genus$Indicator) #155 for indicator column
write.csv(D28_PRRSV_N.genus, file = "D28_PRRSV_N.genus.csv")

D28.genus <- filter(fobar2.gather, grepl(paste("D28"), Day))
summary(D28.genus$Treatment) #Tonsil and Nasal
#   BB Control     IAV   PRRSV 
# 3620    2353    2353    3258 






#Day 36 Control nasal total Genus
D36_Control_N.genus <- filter(fobar2.gather, grepl(paste("D36_Control_N"), All))
D36_Control_N.genus
D36_Control_N.genus$Indicator<-"Zero"
D36_Control_N.genus$Indicator[D36_Control_N.genus$value>0]<-D36_Control_N.genus$Genus[D36_Control_N.genus$value>0]
unique(D36_Control_N.genus$Indicator) #153 for indicator column
write.csv(D36_Control_N.genus, file = "D36_Control_N.genus.csv")

#Day 36 BB nasal total Genus
D36_BB_N.genus <- filter(fobar2.gather, grepl(paste("D36_BB_N"), All))
D36_BB_N.genus
D36_BB_N.genus$Indicator<-"Zero"
D36_BB_N.genus$Indicator[D36_BB_N.genus$value>0]<-D36_BB_N.genus$Genus[D36_BB_N.genus$value>0]
unique(D36_BB_N.genus$Indicator) #158 for indicator column
write.csv(D36_BB_N.genus, file = "D36_BB_N.genus.csv")

#Day 36 IAV nasal total Genus
D36_IAV_N.genus <- filter(fobar2.gather, grepl(paste("D36_IAV_N"), All))
D36_IAV_N.genus
D36_IAV_N.genus$Indicator<-"Zero"
D36_IAV_N.genus$Indicator[D36_IAV_N.genus$value>0]<-D36_IAV_N.genus$Genus[D36_IAV_N.genus$value>0]
unique(D36_IAV_N.genus$Indicator) #154 for indicator column
write.csv(D36_IAV_N.genus, file = "D36_IAV_N.genus.csv")

#Day 36 PRRSV nasal total Genus
D36_PRRSV_N.genus <- filter(fobar2.gather, grepl(paste("D36_PRRSV_N"), All))
D36_PRRSV_N.genus
D36_PRRSV_N.genus$Indicator<-"Zero"
D36_PRRSV_N.genus$Indicator[D36_PRRSV_N.genus$value>0]<-D36_PRRSV_N.genus$Genus[D36_PRRSV_N.genus$value>0]
unique(D36_PRRSV_N.genus$Indicator) #163 for indicator column
write.csv(D36_PRRSV_N.genus, file = "D36_PRRSV_N.genus.csv")

D36.genus <- filter(fobar2.gather, grepl(paste("D36"), Day))
summary(D36.genus$Treatment) #Tonsil and Nasal
#   BB Control     IAV   PRRSV 
# 3620    3439    3258    3620 






#Day 42 Control nasal total Genus
D42_Control_N.genus <- filter(fobar2.gather, grepl(paste("D42_Control_N"), All))
D42_Control_N.genus
D42_Control_N.genus$Indicator<-"Zero"
D42_Control_N.genus$Indicator[D42_Control_N.genus$value>0]<-D42_Control_N.genus$Genus[D42_Control_N.genus$value>0]
unique(D42_Control_N.genus$Indicator) #145 for indicator column
write.csv(D42_Control_N.genus, file = "D42_Control_N.genus.csv")

#Day 42 BB nasal total Genus
D42_BB_N.genus <- filter(fobar2.gather, grepl(paste("D42_BB_N"), All))
D42_BB_N.genus
D42_BB_N.genus$Indicator<-"Zero"
D42_BB_N.genus$Indicator[D42_BB_N.genus$value>0]<-D42_BB_N.genus$Genus[D42_BB_N.genus$value>0]
unique(D42_BB_N.genus$Indicator) #163 for indicator column
write.csv(D42_BB_N.genus, file = "D42_BB_N.genus.csv")

#Day 42 IAV nasal total Genus
D42_IAV_N.genus <- filter(fobar2.gather, grepl(paste("D42_IAV_N"), All))
D42_IAV_N.genus
D42_IAV_N.genus$Indicator<-"Zero"
D42_IAV_N.genus$Indicator[D42_IAV_N.genus$value>0]<-D42_IAV_N.genus$Genus[D42_IAV_N.genus$value>0]
unique(D42_IAV_N.genus$Indicator) #156 for indicator column
write.csv(D42_IAV_N.genus, file = "D42_IAV_N.genus.csv")

#Day 42 PRRSV nasal total Genus
D42_PRRSV_N.genus <- filter(fobar2.gather, grepl(paste("D42_PRRSV_N"), All))
D42_PRRSV_N.genus
D42_PRRSV_N.genus$Indicator<-"Zero"
D42_PRRSV_N.genus$Indicator[D42_PRRSV_N.genus$value>0]<-D42_PRRSV_N.genus$Genus[D42_PRRSV_N.genus$value>0]
unique(D42_PRRSV_N.genus$Indicator) #151 for indicator column
write.csv(D42_PRRSV_N.genus, file = "D42_PRRSV_N.genus.csv")

D42.genus <- filter(fobar2.gather, grepl(paste("D42"), Day))
summary(D42.genus$Treatment) #Tonsil and Nasal
#   BB Control     IAV   PRRSV 
# 3620    3620    3258    3620 