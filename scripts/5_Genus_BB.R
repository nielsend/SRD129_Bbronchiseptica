############################################################
#SRD129 16S - Total Nasal Genera


#NOTES: 
#This code determines the total number of nasal genera found in each treatment group per day and plots as bar graphs
#It also assesses total number of unique families in each day+treatment group 

#Clear workspace and load necessary packages
rm(list=ls())

#Files needed:
#Mothur subsample.shared file: SRD129BB.outsingletons.abund.opti_mcc.0.03.shared
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

otu <- import_mothur(mothur_shared_file = './SRD129BB.outsingletons.abund.opti_mcc.shared')
sample_sums(otu)
taxo <- import_mothur(mothur_constaxonomy_file = './SRD129BB.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
# shared <- read.table('./SRD129BB.outsingletons.abund.opti_mcc.0.03.cons.taxonomy', header = TRUE)
meta <- read.table('./metadata.csv', header = TRUE, sep = ",")
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
sample_sums(SRD129) #sum of all OTUs for each sample
SRD129 <- prune_samples(sample_sums(SRD129) > 7394, SRD129)  #This removes samples that have fewer than 2000 sequences associated with them.
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

fobar.gather <- fobar %>% gather(Genus, value, -(group:Treatment))  #converts to long form dataframe, this is handy for using ggplot faceting functions, check out tidyverse tutorials.
#Create new columns Genus, value; add columns group to All before Genus and value
head(fobar.gather)

#modified ggplot
ggplot(data=fobar.gather, aes(x=Treatment, y=value, fill=Genus)) +
  geom_bar(stat = 'identity') +
  facet_grid(Treatment~Day) + ylab('Percent of total genus in BB and Control groups')

#check if there is an extra "group" column. If so, run the next command and modify which column to remove
which(colnames(phyla_tab2)=="group") #result says column 438
phyla_tab3 <- phyla_tab2[,-438] ##drop the 438th column because it's an extra group (sample name) column
# phyla_tab4 <- phyla_tab3[,colSums(phyla_tab3)>0.1] #keep the columns that have greater than 0.1 value
phyla_tab4 <- phyla_tab3[,colSums(phyla_tab3)>0] #keep the columns that have greater than 0 value

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
unique(fobar2.gather$Genus) #Total number of unique genera: 437
#OR
fobar2.gather %>% summarise_each(funs(n_distinct)) #437 total unique genus

#Calculate percent abundance of each genus per "All" type
fobar2.gather <- fobar2.gather %>% group_by(All) %>% mutate(value2=(value/(length(All)/437))*100)
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
# #For example:
# D28_Control.genus <- subset(fobar2.gather, All=="D28_Control")
# write.csv(D28_Control.genus, file="D28_Control.genus.csv")  
# D28_BB.genus <- subset(fobar2.gather, All=="D28_BB")
# write.csv(D28_BB.genus, file="D28_BB.genus.csv")

#Calculate the total % percent abundance of each genera on each sample (I used JMP to do this) 
#In JMP: copy and paste values from Genus and value columns
#Tables > Summary > input columns for Statistics (% of Total of value column) and Group (Genus column) > OK
#Copy results to a spreadsheet editor (see D0_Control.genus.xlsx for an example)
#Since we are only interested in genera that are above 1% abundance, 
#calculate total percentage of all other genera that are less than 1% in the spreadsheet and label as "Other".
#Create a new spreadsheet and label as "Nasal_genus.csv"
#Create the following columns: Day, Treatment, Genus, Percent Abundance
#Copy the list of genera and their percent abundances from each of the individual Excel files to "Nasal_genus.csv" spreadsheet.
#Fill in the other columns manually (Day, Treatment). 
#You should have a file similar to SRD129_BBControl_NasalGenusPercentAbundance.csv. Continue to the next step.

#Subset each "All" group from fobar2.gather and save as individual csv files.
#Calculate the total % percent abundance of each genera on each sample
#calculate total percentage of all other genera that are less than 1% in the spreadsheet and label as "Other".
#Create a new spreadsheet and label as "nasalgen.csv".

#day 0
D0_Control.genus <- subset(fobar2.gather, All=="D0_Control")
write.csv(D0_Control.genus, file="D0_Control.genus.csv")
D0_Control_sum <- data.frame(aggregate(D0_Control.genus$value, by=list(Category=D0_Control.genus$Genus), FUN=mean))
D0_Control_sum$avg <- (as.numeric(D0_Control_sum$x*100))
D0_Control_sum$Genus[D0_Control_sum$avg < 1] <- "Other"
D0_Control_sum$Genus[D0_Control_sum$avg >= 1] <- D0_Control_sum$Category[D0_Control_sum$avg >= 1]
D0_Control_sum$Treatment <- c("Control") 
D0_Control_sum$Day <- c("D0")
sum(D0_Control_sum$avg)

remove(nasalgen)
nasalgen <- D0_Control_sum

D0_BB.genus <- subset(fobar2.gather, All=="D0_BB")
write.csv(D0_BB.genus, file="D0_BB.genus.csv")
D0_BB_sum <- data.frame(aggregate(D0_BB.genus$value, by=list(Category=D0_BB.genus$Genus), FUN=mean))
D0_BB_sum$avg <- (as.numeric(D0_BB_sum$x*100))
D0_BB_sum$Genus[D0_BB_sum$avg < 1] <- "Other"
D0_BB_sum$Genus[D0_BB_sum$avg >= 1] <- D0_BB_sum$Category[D0_BB_sum$avg >= 1]
D0_BB_sum$Treatment <- c("BB") 
D0_BB_sum$Day <- c("D0")
sum(D0_BB_sum$`x`)
sum(D0_BB_sum$avg)

nasalgen <- rbind(nasalgen, D0_BB_sum)

#day 1
D1_Control.genus <- subset(fobar2.gather, All=="D1_Control")
write.csv(D1_Control.genus, file="D1_Control.genus.csv")
D1_Control_sum <- data.frame(aggregate(D1_Control.genus$value, by=list(Category=D1_Control.genus$Genus), FUN=mean))
D1_Control_sum$avg <- (as.numeric(D1_Control_sum$x*100))
D1_Control_sum$Genus[D1_Control_sum$avg < 1] <- "Other"
D1_Control_sum$Genus[D1_Control_sum$avg >= 1] <- D1_Control_sum$Category[D1_Control_sum$avg >= 1]
D1_Control_sum$Treatment <- c("Control") 
D1_Control_sum$Day <- c("D1")
sum(D1_Control_sum$avg)
nasalgen <- rbind(nasalgen, D1_Control_sum)


D1_BB.genus <- subset(fobar2.gather, All=="D1_BB")
write.csv(D1_BB.genus, file="D1_BB.genus.csv")
D1_BB_sum <- data.frame(aggregate(D1_BB.genus$value, by=list(Category=D1_BB.genus$Genus), FUN=mean))
D1_BB_sum$avg <- (as.numeric(D1_BB_sum$x*100))
D1_BB_sum$Genus[D1_BB_sum$avg < 1] <- "Other"
D1_BB_sum$Genus[D1_BB_sum$avg >= 1] <- D1_BB_sum$Category[D1_BB_sum$avg >= 1]
D1_BB_sum$Treatment <- c("BB") 
D1_BB_sum$Day <- c("D1")
sum(D1_BB_sum$`x`)
sum(D1_BB_sum$avg)

nasalgen <- rbind(nasalgen, D1_BB_sum)





#day 3
D3_BB.genus <- subset(fobar2.gather, All=="D3_BB")
write.csv(D3_BB.genus, file="D3_BB.genus.csv")
D3_BB_sum <- data.frame(aggregate(D3_BB.genus$value, by=list(Category=D3_BB.genus$Genus), FUN=mean))
D3_BB_sum$avg <- (as.numeric(D3_BB_sum$x*100))
D3_BB_sum$Genus[D3_BB_sum$avg < 1] <- "Other"
D3_BB_sum$Genus[D3_BB_sum$avg >= 1] <- D3_BB_sum$Category[D3_BB_sum$avg >= 1]
D3_BB_sum$Treatment <- c("BB") 
D3_BB_sum$Day <- c("D3")
sum(D3_BB_sum$`x`)
sum(D3_BB_sum$avg)


nasalgen <- rbind(nasalgen, D3_BB_sum)


D3_Control.genus <- subset(fobar2.gather, All=="D3_Control")
write.csv(D3_Control.genus, file="D3_Control.genus.csv")
D3_Control_sum <- data.frame(aggregate(D3_Control.genus$value, by=list(Category=D3_Control.genus$Genus), FUN=mean))
D3_Control_sum$avg <- (as.numeric(D3_Control_sum$x*100))
D3_Control_sum$Genus[D3_Control_sum$avg < 1] <- "Other"
D3_Control_sum$Genus[D3_Control_sum$avg >= 1] <- D3_Control_sum$Category[D3_Control_sum$avg >= 1]
D3_Control_sum$Treatment <- c("Control") 
D3_Control_sum$Day <- c("D3")
sum(D3_Control_sum$avg)

nasalgen <- rbind(nasalgen, D3_Control_sum)



#day 7
D7_BB.genus <- subset(fobar2.gather, All=="D7_BB")
write.csv(D7_BB.genus, file="D7_BB.genus.csv")
D7_BB_sum <- data.frame(aggregate(D7_BB.genus$value, by=list(Category=D7_BB.genus$Genus), FUN=mean))
D7_BB_sum$avg <- (as.numeric(D7_BB_sum$x*100))
D7_BB_sum$Genus[D7_BB_sum$avg < 1] <- "Other"
D7_BB_sum$Genus[D7_BB_sum$avg >= 1] <- D7_BB_sum$Category[D7_BB_sum$avg >= 1]
D7_BB_sum$Treatment <- c("BB") 
D7_BB_sum$Day <- c("D7")
sum(D7_BB_sum$`x`)
sum(D7_BB_sum$avg)

nasalgen <- rbind(nasalgen, D7_BB_sum)


D7_Control.genus <- subset(fobar2.gather, All=="D7_Control")
write.csv(D7_Control.genus, file="D7_Control.genus.csv")
D7_Control_sum <- data.frame(aggregate(D7_Control.genus$value, by=list(Category=D7_Control.genus$Genus), FUN=mean))
D7_Control_sum$avg <- (as.numeric(D7_Control_sum$x*100))
D7_Control_sum$Genus[D7_Control_sum$avg < 1] <- "Other"
D7_Control_sum$Genus[D7_Control_sum$avg >= 1] <- D7_Control_sum$Category[D7_Control_sum$avg >= 1]
D7_Control_sum$Treatment <- c("Control") 
D7_Control_sum$Day <- c("D7")
sum(D7_Control_sum$avg)
nasalgen <- rbind(nasalgen, D7_Control_sum)



#day 10
D10_BB.genus <- subset(fobar2.gather, All=="D10_BB")
write.csv(D10_BB.genus, file="D10_BB.genus.csv")
D10_BB_sum <- data.frame(aggregate(D10_BB.genus$value, by=list(Category=D10_BB.genus$Genus), FUN=mean))
D10_BB_sum$avg <- (as.numeric(D10_BB_sum$x*100))
D10_BB_sum$Genus[D10_BB_sum$avg < 1] <- "Other"
D10_BB_sum$Genus[D10_BB_sum$avg >= 1] <- D10_BB_sum$Category[D10_BB_sum$avg >= 1]
D10_BB_sum$Treatment <- c("BB") 
D10_BB_sum$Day <- c("D10")
sum(D10_BB_sum$`x`)
sum(D10_BB_sum$avg)


nasalgen <- rbind(nasalgen, D10_BB_sum)


D10_Control.genus <- subset(fobar2.gather, All=="D10_Control")
write.csv(D10_Control.genus, file="D10_Control.genus.csv")
D10_Control_sum <- data.frame(aggregate(D10_Control.genus$value, by=list(Category=D10_Control.genus$Genus), FUN=mean))
D10_Control_sum$avg <- (as.numeric(D10_Control_sum$x*100))
D10_Control_sum$Genus[D10_Control_sum$avg < 1] <- "Other"
D10_Control_sum$Genus[D10_Control_sum$avg >= 1] <- D10_Control_sum$Category[D10_Control_sum$avg >= 1]
D10_Control_sum$Treatment <- c("Control") 
D10_Control_sum$Day <- c("D10")
sum(D10_Control_sum$avg)

nasalgen <- rbind(nasalgen, D10_Control_sum)



#day 14
D14_BB.genus <- subset(fobar2.gather, All=="D14_BB")
write.csv(D14_BB.genus, file="D14_BB.genus.csv")
D14_BB_sum <- data.frame(aggregate(D14_BB.genus$value, by=list(Category=D14_BB.genus$Genus), FUN=mean))
D14_BB_sum$avg <- (as.numeric(D14_BB_sum$x*100))
D14_BB_sum$Genus[D14_BB_sum$avg < 1] <- "Other"
D14_BB_sum$Genus[D14_BB_sum$avg >= 1] <- D14_BB_sum$Category[D14_BB_sum$avg >= 1]
D14_BB_sum$Treatment <- c("BB") 
D14_BB_sum$Day <- c("D14")
sum(D14_BB_sum$`x`)
sum(D14_BB_sum$avg)

nasalgen <- rbind(nasalgen, D14_BB_sum)


D14_Control.genus <- subset(fobar2.gather, All=="D14_Control")
write.csv(D14_Control.genus, file="D14_Control.genus.csv")
D14_Control_sum <- data.frame(aggregate(D14_Control.genus$value, by=list(Category=D14_Control.genus$Genus), FUN=mean))
D14_Control_sum$avg <- (as.numeric(D14_Control_sum$x*100))
D14_Control_sum$Genus[D14_Control_sum$avg < 1] <- "Other"
D14_Control_sum$Genus[D14_Control_sum$avg >= 1] <- D14_Control_sum$Category[D14_Control_sum$avg >= 1]
D14_Control_sum$Treatment <- c("Control") 
D14_Control_sum$Day <- c("D14")
sum(D14_Control_sum$avg)
nasalgen <- rbind(nasalgen, D14_Control_sum)



#day 21
D21_BB.genus <- subset(fobar2.gather, All=="D21_BB")
write.csv(D21_BB.genus, file="D21_BB.genus.csv")
D21_BB_sum <- data.frame(aggregate(D21_BB.genus$value, by=list(Category=D21_BB.genus$Genus), FUN=mean))
D21_BB_sum$avg <- (as.numeric(D21_BB_sum$x*100))
D21_BB_sum$Genus[D21_BB_sum$avg < 1] <- "Other"
D21_BB_sum$Genus[D21_BB_sum$avg >= 1] <- D21_BB_sum$Category[D21_BB_sum$avg >= 1]
D21_BB_sum$Treatment <- c("BB") 
D21_BB_sum$Day <- c("D21")
sum(D21_BB_sum$`x`)
sum(D21_BB_sum$avg)

nasalgen <- rbind(nasalgen, D21_BB_sum)


D21_Control.genus <- subset(fobar2.gather, All=="D21_Control")
write.csv(D21_Control.genus, file="D21_Control.genus.csv")
D21_Control_sum <- data.frame(aggregate(D21_Control.genus$value, by=list(Category=D21_Control.genus$Genus), FUN=mean))
D21_Control_sum$avg <- (as.numeric(D21_Control_sum$x*100))
D21_Control_sum$Genus[D21_Control_sum$avg < 1] <- "Other"
D21_Control_sum$Genus[D21_Control_sum$avg >= 1] <- D21_Control_sum$Category[D21_Control_sum$avg >= 1]
D21_Control_sum$Treatment <- c("Control") 
D21_Control_sum$Day <- c("D21")
sum(D21_Control_sum$avg)
nasalgen <- rbind(nasalgen, D21_Control_sum)




#day 36
D36_BB.genus <- subset(fobar2.gather, All=="D36_BB")
write.csv(D36_BB.genus, file="D36_BB.genus.csv")
D36_BB_sum <- data.frame(aggregate(D36_BB.genus$value, by=list(Category=D36_BB.genus$Genus), FUN=mean))
D36_BB_sum$avg <- (as.numeric(D36_BB_sum$x*100))
D36_BB_sum$Genus[D36_BB_sum$avg < 1] <- "Other"
D36_BB_sum$Genus[D36_BB_sum$avg >= 1] <- D36_BB_sum$Category[D36_BB_sum$avg >= 1]
D36_BB_sum$Treatment <- c("BB") 
D36_BB_sum$Day <- c("D36")
sum(D36_BB_sum$`x`)
sum(D36_BB_sum$avg)

nasalgen <- rbind(nasalgen, D36_BB_sum)



D36_Control.genus <- subset(fobar2.gather, All=="D36_Control")
write.csv(D36_Control.genus, file="D36_Control.genus.csv")
D36_Control_sum <- data.frame(aggregate(D36_Control.genus$value, by=list(Category=D36_Control.genus$Genus), FUN=mean))
D36_Control_sum$avg <- (as.numeric(D36_Control_sum$x*100))
D36_Control_sum$Genus[D36_Control_sum$avg < 1] <- "Other"
D36_Control_sum$Genus[D36_Control_sum$avg >= 1] <- D36_Control_sum$Category[D36_Control_sum$avg >= 1]
D36_Control_sum$Treatment <- c("Control") 
D36_Control_sum$Day <- c("D36")
sum(D36_Control_sum$avg)

nasalgen <- rbind(nasalgen, D36_Control_sum)


#day 42
D42_BB.genus <- subset(fobar2.gather, All=="D42_BB")
write.csv(D42_BB.genus, file="D42_BB.genus.csv")
D42_BB_sum <- data.frame(aggregate(D42_BB.genus$value, by=list(Category=D42_BB.genus$Genus), FUN=mean))
D42_BB_sum$avg <- (as.numeric(D42_BB_sum$x*100))
D42_BB_sum$Genus[D42_BB_sum$avg < 1] <- "Other"
D42_BB_sum$Genus[D42_BB_sum$avg >= 1] <- D42_BB_sum$Category[D42_BB_sum$avg >= 1]
D42_BB_sum$Treatment <- c("BB") 
D42_BB_sum$Day <- c("D42")
sum(D42_BB_sum$`x`)
sum(D42_BB_sum$avg)

nasalgen <- rbind(nasalgen, D42_BB_sum)



D42_Control.genus <- subset(fobar2.gather, All=="D42_Control")
write.csv(D42_Control.genus, file="D42_Control.genus.csv")
D42_Control_sum <- data.frame(aggregate(D42_Control.genus$value, by=list(Category=D42_Control.genus$Genus), FUN=mean))
D42_Control_sum$avg <- (as.numeric(D42_Control_sum$x*100))
D42_Control_sum$Genus[D42_Control_sum$avg < 1] <- "Other"
D42_Control_sum$Genus[D42_Control_sum$avg >= 1] <- D42_Control_sum$Category[D42_Control_sum$avg >= 1]
D42_Control_sum$Treatment <- c("Control") 
D42_Control_sum$Day <- c("D42")
sum(D42_Control_sum$avg)

nasalgen <- rbind(nasalgen, D42_Control_sum)
write.csv(nasalgen, "./nasalgen.csv")


##### Plot nasal genera above 1% abundance ######
#nasalgen <- read.table('./data/SRD129_BBControl_NasalGenusPercentAbundance.csv', header = TRUE, sep = ",")
head(nasalgen)
unique(nasalgen$Treatment) #Control BB
unique(nasalgen$Day) #D0  D1  D10 D14 D21 D28 D3  D36 D42 D7
nasalgen$Day = factor(nasalgen$Day, levels = c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D36", "D42"))
levels(nasalgen$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21"  "D36" "D42"
nasalgenplot <- ggplot(data=nasalgen, aes(x=Treatment, y=avg, fill=Genus)) +
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
  scale_fill_manual(values = mokoleColor, aesthetics = "fill")
nasalgenplot <- nasalgenplot + guides(fill= guide_legend(ncol = 2))
nasalgenplot <- nasalgenplot + theme(legend.text = element_text(face = 'italic')) +
  scale_fill_manual(values = mokoleColor, aesthetics = "fill")
nasalgenplot
unique(nasalgen$Genus)

#Save 'nasalgenplot' as a .tiff for publication, 500dpi
ggsave("nasalgenplot.tiff", plot=nasalgenplot, width = 15, height = 6, dpi = 500, units =c("in"))


 mokoleColor <- c("#a9a9a9", "#2f4f4f", "#8b4513", "#2e8b57",
                 "#191970", "#006400", "#8b0000", "#808000", "#008b8b", "#4682b4",
                  "#9acd32", "#00008b", "#daa520", "#7f007f", "#66cdaa", "#9932cc",
                  "#ff0000", "#ff8c00", "#c71585", "#00ff7f",
                 "#e9967a", "#dc143c", "#00ffff", "#00bfff", "#9370db",
                "#0000ff", "#adff2f", "#ff7f50", "#ff00ff", "#1e90ff",
                 "#db7093", "#FF5733", "#ffff54", "#dda0dd", "#90ee90", 
                "#afeeee", "#ee82ee",  "#5C4033")

# mokoleColor <- c("#a9a9a9", "#2f4f4f", "#8b4513", "#2e8b57",
#                  "#191970", "#006400", "#8b0000", "#808000", "#008b8b", "#4682b4",
#                  "#9acd32", "#00008b", "#daa520", "#7f007f", "#66cdaa", "#9932cc",
#                  "#ff0000", "#ff8c00", "#c71585", "#00ff00", "#00ff7f",
#                  "#e9967a", "#dc143c", "#00ffff", "#00bfff", "#9370db",
#                  "#0000ff", "#adff2f", "#ff7f50", "#ff00ff", "#1e90ff",
#                  "#db7093", "#f0e68c", "#ffff54", "#dda0dd", "#90ee90", 
#                  "#afeeee", "#ee82ee", "#ffe4c4", "#ffc0cb")

 
nasalgen$Genus <- factor(nasalgen$Genus, levels = c("Acinetobacter", "Actinobacillus", "Aerococcus", "Alloprevotella", 
                                                    "Bergeyella", "Blautia", "Bordetella","Carnobacteriaceae_unclassified", 
                                                   "Chryseobacterium", "Clostridium_sensu_stricto_1", "Empedobacter",
                                                   "Enhydrobacter", "Escherichia-Shigella", "Kurthia", "Lachnospiraceae_unclassified", 
                                                   "Lachnospiraceae_XPB1014_group", "Lactobacillus", "Lactococcus", "Megasphaera",
                                                   "Moraxella", "Neisseriaceae_unclassified", "Pasteurellaceae_unclassified",
                                                   "Phascolarctobacterium", "Prevotella_9", "Prevotellaceae_NK3B31_group", 
                                                   "Prevotellaceae_UCG-003", "Prevotellaceae_unclassified", 
                                                   "Rikenellaceae_RC9_gut_group", "Rothia", "Ruminococcaceae_UCG-005", 
                                                   "Staphylococcus","Streptococcus", "Terrisporobacter", 
                                                   "Treponema_2", "uncultured", "Weeksellaceae_unclassified", "Weissella", "Other"))
nasalgen %>% filter(Genus=="Bordetella") 
View(taxonomy)
taxonomy2 <- taxonomy %>% separate(Taxonomy, sep=";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
taxonomy2 <- taxonomy2 %>% filter(Genus=="Bordetella")

BordetellaGroupDays <- fobar2.gather %>% filter(Genus=="Bordetella" & Day != "D0" & Treatment=="BB")
write.csv(BordetellaGroupDays, "./BordetellaGroupDays.csv")

BGroupDays <- phyla_tab$Bordetella
BGroupDays<- as.data.frame(BGroupDays)
BGroupDays$group <- row.names(phyla_tab)
BGroupMerge <- merge(BGroupDays, meta, by="group")
BGroupMerge <- BGroupMerge %>% filter(Day != "D0" & Treatment=="BB")

write.csv(BGroupMerge, "./BGroupMerge.csv")

colnames(phyla_tab)



