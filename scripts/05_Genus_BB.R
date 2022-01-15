############################################################
#SRD129 16S - Total Genera
#By Mou KT

#NOTES: 
#This code determines the total number of nasal genera found in each treatment group per day and plots as bar graphs
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

#write csv and save in ./data
write.csv(fobar2.gather, file="SRD129_BBControl_NasalGenusPercentAbundance.csv")

#### Subset each "All" group from fobar2.gather and save as individual csv files ###
#For example:
D28_Control.genus <- subset(fobar2.gather, All=="D28_Control")
write.csv(D28_Control.genus, file="D28_Control.genus.csv")  
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
#You should have a file similar to SRD129_BBControl_NasalGenusPercentAbundance.csv. Continue to the next step.

##### Plot nasal genera above 1% abundance ######
nasalgen <- read.table('./data/SRD129_BBControl_NasalGenusPercentAbundance.csv', header = TRUE, sep = ",")
head(nasalgen)
unique(nasalgen$Treatment) #Control BB
unique(nasalgen$Day) #D0  D1  D10 D14 D21 D28 D3  D36 D42 D7
nasalgen$Day = factor(nasalgen$Day, levels = c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
levels(nasalgen$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
nasalgenplot <- ggplot(data=nasalgen, aes(x=Treatment, y=Percent.abundance, fill=Genus)) +
    geom_bar(stat = 'identity') +
    facet_grid(~Day) + ylab('Relative abundance (%) at nasal site') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(size=12, angle=45, hjust=1),
          axis.title.x = element_blank(),
          strip.text = element_text(size= 15),
          axis.text.y = element_text(size=14), 
          axis.title.y = element_text(size=14), 
          legend.text=element_text(size=14), 
          legend.title=element_text(size=14)) +
    scale_fill_igv(name = "Genus")
nasalgenplot <- nasalgenplot + guides(fill= guide_legend(ncol = 2))
nasalgenplot <- nasalgenplot + theme(legend.text = element_text(face = 'italic'))
nasalgenplot