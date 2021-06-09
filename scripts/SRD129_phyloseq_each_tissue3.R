#########################################
#SRD129 16S phyloseq prep, adonis, PCoA each tissue- BB and control
#Kathy Mou

#NOTES: 
#This code creates phyloseq objects ("All" or no "All" column) for beta diversity statistics 
#of nasal samples
## a couple adonis functions and results that you should include in paper (treatment,
#day, day&treatment effects on variation)
#PCoA plots (day, treatment, tissue if applicable, treatment x day) to see general trends for each tissue

#Clear workspace and load necessary packages
rm(list=ls())

#Set working directory on desktop
#Mac
setwd("~/Desktop/Microbiome/Projects/SRD129/SRD129_2000singletons")
setwd("~/Desktop/SRD129/SRD129_2000singletons")

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)

#Load image file
load("SRD129_phyloseq_each_tissue_nobad3.RData")


#Save image file
save.image(file="SRD129_phyloseq_each_tissue_nobad3.RData")
#########################################

#Read files
otu <- read.csv("SRD129abundsingleton2000OTUtable.csv", row.names=1) #set column 1 as row names
meta <-read.csv("SRD129metadata06292018.csv")

dim(otu) #2196 909
head(otu[,1:5])
head(otu[,900:909])

dim(meta) #904 7
head(meta[,1:7])

#Remove taxonomy
tax <- otu[,(908:909)] #removed column 908 taxonomy and 909 to include the row names
head(tax)
colnames(tax)[1] <- "delete" #renamed column 1 (formerly 908) as "delete" which will later be deleted
head(tax)

#OTU count data only
otu <- otu[,-909] #remove column 909 taxonomy to have only otu data
head(otu[,900:908]) 
dim(otu) #2196 rows 908 columns

#Transpose to match metadata format
otu.trans <- t(otu) #now rownames are sample names, columns are OTUs
head(otu.trans[,1:5])
head(meta) #row names are numbered, but want sample names as row names
meta <- meta[,-1]
head(meta)

class(meta) #dataframe
class(otu) #dataframe

#Merge otu and meta data frames
otu.meta <- merge(meta, otu.trans, by.x=1, by.y=0) 
#x=1 means match via 1st column from meta; y=0 means match via rownames from otu.trans
head(otu.meta[,1:10])
class(otu.meta) #dataframe
dim(otu.meta) #904 2202

#With "All" column
write.csv(otu.meta, file="SRD129abundsingleton2000.otu.meta_all.csv")

#Subset tissues with "All" column
nw <- subset(otu.meta, Tissue=="N")
head(nw[,1:10])
dim(nw) #472 rows 2202 columns

write.csv(nw, file="SRD129abundsingleton2000.otu.meta_all.nasalonly.csv")

tt <- subset(otu.meta, Tissue=="T")
head(tt[,1:10])
tail(tt[,1:10])
dim(tt) #431 rows 2202 columns

write.csv(tt, file="SRD129abundsingleton2000.otu.meta_all.tonsilonly.csv")

#Searching for the sample that isn't just an N or T
otu.meta.subset <- otu.meta[grep("NT", otu.meta$Tissue), ]
head(otu.meta.subset[,1:5])
#Pig 132, day 1, NT
#Did not include this in with Nasal or Tonsil only analysis, but will make a separate file
write.csv(otu.meta.subset, file="SRD129abundsingleton2000.otu.meta_all.P132D1NT.csv")

#--------------TONSIL PHYLOSEQ with "All" column--------------#

#Create tonsil (with "All" column) phyloseq object
library("phyloseq")
library("tidyr")

#meta only
head(tt[,1:10])
meta.tt <- tt[,1:6] #take columns 1-6 (Sample to All) to make meta.tt
head(meta.tt)
row.names(meta.tt) <- meta.tt[,1] #make column 1 be rownames for meta.tt
head(meta.tt)
meta.tt <- meta.tt[,-1] #remove the extra "Sample" column
head(meta.tt)
dim(meta.tt) #431 5

#Create SAM metada table phyloseq object
SAMtt = sample_data(meta.tt, errorIfNULL = TRUE)
head(SAMtt)
dim(SAMtt) #431 5

#OTU only
head(tt[,1:10])
dim(tt) #431 2202
row.names(tt) <- tt[,1] #make column 1 be rownames for tt
head(tt[,1:10])
tt <- tt[,-1] #remove the extra "Sample" column
head(tt[,1:10])
dim(tt) #431 2201
head(tt[,2199:2201])
otu.tt <- tt[,c(6:2201)] #select Sample and Otu columns to create otu.tt dataframe
head(otu.tt[,1:10])
dim(otu.tt) #431 2196
otu.tt.trans <- t(otu.tt) #transpose otu.tt to have Otu as rownames, sample names as column names
head(otu.tt.trans[,1:10])
dim(otu.tt.trans) #2196 431

#Merge tax back into otu for correct format and taxons
head(tax)
otu.tax.tt <- merge(otu.tt.trans, tax, by=0) #merge by rownames aka Otu rownames
dim(otu.tax.tt) #2196 434
head(otu.tax.tt[,430:434])
head(otu.tax.tt[,1:5])
row.names(otu.tax.tt) <- otu.tax.tt[,1] #set first row as rownames
head(otu.tax.tt[,1:5])
otu.tax.tt <- otu.tax.tt[,-1] #remove first row, extraneious Otu column
head(otu.tax.tt[,1:5])

#Split again
dim(otu.tax.tt)#2196 433
head(otu.tax.tt[,430:433])
otu.notax.tt <- otu.tax.tt[,1:431] #take rows 1-431 to make new dataframe otu.notax.tt2 (36 is delete column, 37 is taxonomy column)
head(otu.notax.tt[,1:5])
head(otu.notax.tt[,425:431])
dim(otu.notax.tt) #2196 431
class(otu.notax.tt) #dataframe
otu.notax.tt <- as.matrix(otu.notax.tt) #turn otu.notax.tt2 into a matrix class
class(otu.notax.tt) #matrix

#Create OTU table phyloseq object
OTUtt = otu_table(otu.notax.tt, taxa_are_rows = TRUE)
head(OTUtt)
dim(OTUtt) #2196 431

head(otu.tax.tt[,425:433])
tax.levels.tt <- separate(data = otu.tax.tt, 
                       col = Taxonomy, 
                       into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#Separate Taxonomy column into 7 separate columns labeled Kingdom to Species
head(tax.levels.tt) #notice that Species column is blank
dim(tax.levels.tt) #2196 439
head(tax.levels.tt[,430:439])
tax.only.tt <- tax.levels.tt[,433:439] #keep only taxonomy columns Kingdom to Genus
head(tax.only.tt)
dim(tax.only.tt) #2196 6
class(tax.only.tt) #dataframe
tax.m.tt <- as.matrix(tax.only.tt)
class(tax.m.tt) #matrix

#Create TAX taxonomy table phyloseq object
TAXtt = tax_table(tax.m.tt)
head(TAXtt)
dim(TAXtt) #2196 7
TAXtt <- TAXtt[,-7]
head(TAXtt)
dim(TAXtt) #2196 6

#Create phyloseq object containing taxonomy, metadata, and otu table
phyloseq.tt <- phyloseq(OTUtt, SAMtt, TAXtt)
phyloseq.tt #view phyloseq object
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 2196 taxa and 431 samples ]
#sample_data() Sample Data:       [ 431 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 2196 taxa by 6 taxonomic ranks ]

save(phyloseq.tt, file="SRD129.phyloseq.tt.RData")



##-------NASAL PHYLOSEQ with "All" column--------##

#Create nasal phyloseq object
library("phyloseq")

#Nasal
#Split info again
#meta only
head(nw[,1:10])
dim(nw) #472 2202
meta.nw <- nw[,1:6] #take columns 1-6 to make meta.nw
head(meta.nw)
row.names(meta.nw) <- meta.nw[,1] #make column 1 be rownames for meta.nw
head(meta.nw)
meta.nw <- meta.nw[,-1] #remove the extra "Sample" column
head(meta.nw)
dim(meta.nw) #472 5
class(meta.nw)

#Create SAM metada table phyloseq object
SAMnw = sample_data(meta.nw, errorIfNULL = TRUE)
head(SAMnw)
dim(SAMnw) #472 5

#OTU only
head(nw[,1:10])
dim(nw) #472 2202
head(nw[,2199:2202])
otu.nw <- nw[,c(1,7:2202)] #select Sample and Otu columns to create otu.nw dataframe
head(otu.nw[,1:10])
dim(otu.nw) #472 2197
row.names(otu.nw) <- otu.nw[,1] #make first column be rownames for otu.nw
head(otu.nw[,1:10])
otu.nw <- otu.nw[,-1] #remove the first column
head(otu.nw[,1:10])
otu.nw.trans <- t(otu.nw) #transpose otu.nw to have Otu as rownames, sample names as column names
head(otu.nw.trans[,1:10])
dim(otu.nw.trans) #2196 472
head(otu.nw.trans[,469:472])

#Merge tax back into otu for correct format and taxons
head(tax)
otu.tax.nw <- merge(otu.nw.trans, tax, by=0) #merge by rownames aka Otu rownames
dim(otu.tax.nw) #2196 475
head(otu.tax.nw[,470:475])
head(otu.tax.nw[,1:5])
row.names(otu.tax.nw) <- otu.tax.nw[,1] #set first row as rownames
head(otu.tax.nw[,1:5])
otu.tax.nw <- otu.tax.nw[,-1] #remove first row, extraneous Otu column
head(otu.tax.nw[,1:5])
dim(otu.tax.nw) #2196 474

#Split again
dim(otu.tax.nw) #2196 474
head(otu.tax.nw[,470:474])
otu.notax.nw <- otu.tax.nw[,1:472] #take rows 1-472 to make new dataframe otu.notax.nw (473 is delete column, 474 is taxonomy column)
head(otu.notax.nw[,1:5])
head(otu.notax.nw[,468:472])
dim(otu.notax.nw) #2196 472
class(otu.notax.nw) #dataframe
otu.notax.nw <- as.matrix(otu.notax.nw) #turn otu.notax.nw into a matrix class
class(otu.notax.nw) #matrix

#Create OTU table phyloseq object
OTUnw = otu_table(otu.notax.nw, taxa_are_rows = TRUE)
head(OTUnw)
dim(OTUnw) #2196 472
class(OTUnw)

dim(otu.tax.nw) #2196 474
head(otu.tax.nw[,470:474])
tax.levels.nw <- separate(data = otu.tax.nw, 
                        col = Taxonomy, 
                        into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#Separate Taxonomy column into 7 separate columns labeled Kingdom to Species
head(tax.levels.nw) #notice that Species column is blank
dim(tax.levels.nw) #2196 480
head(tax.levels.nw[,473:480])
tax.only.nw <- tax.levels.nw[,474:480] #keep only taxonomy columns Kingdom to Genus
head(tax.only.nw)
dim(tax.only.nw) #2196 7
class(tax.only.nw) #data.frame
tax.m.nw <- as.matrix(tax.only.nw)
class(tax.m.nw) #matrix
head(tax.m.nw)

#Create TAX taxonomy table phyloseq object
TAXnw = tax_table(tax.m.nw)
head(TAXnw)
dim(TAXnw) #2196 7
TAXnw <- TAXnw[,-7]
dim(TAXnw) #2196 6
head(TAXnw)
class(TAXnw)


#Create phyloseq object containing taxonomy, metadata, and otu table
phyloseq.nw <- phyloseq(OTUnw, SAMnw, TAXnw)
phyloseq.nw #view phyloseq object
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 2196 taxa and 472 samples ]
#sample_data() Sample Data:       [ 472 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 2196 taxa by 6 taxonomic ranks ]

save(phyloseq.nw, file="SRD129.phyloseq.nw.RData")

##############################################################################################################################
#Remove bad samples

#ATTEMPT #1 --- DID NOT WORK WHEN I TRIED RUNNING ADONIS FUNCTION
phyloseq.nw.nobad = subset_samples(phyloseq.nw, row.names(sample_data(phyloseq.nw)) != "SRD129_P150_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P151_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P152_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P153_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P154_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P156_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P163_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P164_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P165_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P166_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P167_N_D28" & row.names(sample_data(phyloseq.nw)) != "SRD129_P168_N_D28")
#to remove specific samples, must include row.names(sample_data(phyloseq.nw)) != "<sample name>"
#add an "&" to add additional samples to exclude
phyloseq.nw.nobad
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 2196 taxa and 460 samples ]
#sample_data() Sample Data:       [ 460 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 2196 taxa by 6 taxonomic ranks ]
head(sample_data(phyloseq.nw.nobad))
tail(sample_data(phyloseq.nw.nobad))

#Save the phyloseq object without bad samples as a separate phyloseq object
save(phyloseq.nw.nobad, file="SRD129.phyloseq.nw.nobad.RData")


#ATTEMPT #2 --- ALSO DID NOT WORK ON ADONIS FUNCTION
#Read files
otunobad <- read.csv("SRD129abundsingleton2000OTUtablenobad.csv", row.names=1) #set column 1 as row names
metanobad <-read.csv("SRD129metadata06292018nobad.csv")

dim(otunobad) #2196 897
head(otunobad[,1:5])
head(otunobad[,890:897])

dim(metanobad) #892 7
head(metanobad[,1:7])

#Remove taxonomy
taxnobad <- otunobad[,(896:897)] #removed column 897 taxonomy and 896 to include the row names
head(taxnobad)
colnames(taxnobad)[1] <- "delete" #renamed column 1 (formerly 896) as "delete" which will later be deleted
head(taxnobad)

#OTU count data only
otunobad <- otunobad[,-897] #remove column 897 taxonomy to have only otu data
head(otunobad)
head(otunobad[,890:896]) 
dim(otunobad) #2196 rows 896 columns

#Transpose to match metadata format
otunobad.trans <- t(otunobad) #now rownames are sample names, columns are OTUs
head(otunobad.trans[,1:5])
head(metanobad) #row names are numbered, but want sample names as row names
metanobad <- metanobad[,-1]
dim(metanobad) #892 6
rownames(otunobad.trans) #contains 4 Mock samples

class(metanobad) #dataframe
class(otunobad) #dataframe

#Merge otu and meta data frames
otunobad.meta <- merge(metanobad, otunobad.trans, by.x=1, by.y=0) 
#x=1 means match via 1st column from meta; y=0 means match via rownames from otunobad.trans
head(otunobad.meta[,1:10])
tail(otunobad.meta)
class(otunobad.meta) #dataframe
dim(otunobad.meta) #892 2202

#With "All" column
write.csv(otunobad.meta, file="SRD129abundsingleton2000.otu.meta_allnobad.csv")

#Subset tissues with "All" column
nw.nobadbb <- subset(otunobad.meta, Tissue=="N")
head(nw.nobadbb[,1:10])
dim(nw.nobadbb) #460 rows 2202 columns

#write.csv(nw.nobad, file="SRD129abundsingleton2000.otu.meta_allnobad.nasalonly.csv")

##-------NASAL PHYLOSEQ with "All" column, no bad samples, BB and control only, no DENG12, DNEG6--------##

#Nasal
#Split info again
#meta only
head(nw.nobadbb[,1:10])
dim(nw.nobadbb) #460 2202
colnames(nw.nobadbb)
rownames(nw.nobadbb)
unique(nw.nobadbb$Treatment)
nw.nobadb<- nw.nobadbb[!grepl("PRRSV", nw.nobadbb$Treatment),]
nw.nobadb<- nw.nobadb[!grepl("IAV", nw.nobadb$Treatment),]
unique(nw.nobadb$Treatment) #BB Control
nw.nobadb <- nw.nobadb[!grepl("DNEG12", nw.nobadb$Day),]
nw.nobadb <- nw.nobadb[!grepl("DNEG6", nw.nobadb$Day),]
unique(nw.nobadb$Day) #D0  D10 D14 D21 D28 D3  D36 D42 D7  D1
head(nw.nobadb[,1:10])


meta.nw.nobad <- nw.nobadb[,1:6] #take columns 1-6 to make meta.nw.nobad
head(meta.nw.nobad)
row.names(meta.nw.nobad) <- meta.nw.nobad[,1] #make column 1 be rownames for meta.nw.nobad
head(meta.nw.nobad)
meta.nw.nobad <- meta.nw.nobad[,-1] #remove the extra "Sample" column
head(meta.nw.nobad)
dim(meta.nw.nobad) #191 5
class(meta.nw.nobad) #data.frame

#Create SAM metada table phyloseq object
SAMnw.nobad = sample_data(meta.nw.nobad, errorIfNULL = TRUE)
head(SAMnw.nobad)
dim(SAMnw.nobad) #191 5

#OTU only
head(nw.nobadb[,1:10])
dim(nw.nobadb) #191 2202
head(nw.nobadb[,2199:2202])
otu.nw.nobad <- nw.nobadb[,c(1,7:2202)] #select Sample and Otu columns to create otu.nw.nobad dataframe
head(otu.nw.nobad[,1:10])
dim(otu.nw.nobad) #191 2197
row.names(otu.nw.nobad) <- otu.nw.nobad[,1] #make first column be rownames for otu.nw.nobad
head(otu.nw.nobad[,1:10])
otu.nw.nobad <- otu.nw.nobad[,-1] #remove the first column
head(otu.nw.nobad[,1:10])
otu.nw.nobad.trans <- t(otu.nw.nobad) #transpose otu.nw.nobad to have Otu as rownames, sample names as column names
head(otu.nw.nobad.trans[,1:10])
dim(otu.nw.nobad.trans) #2196 191
head(otu.nw.nobad.trans[,185:191])

#Merge tax back into otu for correct format and taxons
head(taxnobad)
otu.tax.nw.nobad <- merge(otu.nw.nobad.trans, tax, by=0) #merge by rownames aka Otu rownames
dim(otu.tax.nw.nobad) #2196 194
head(otu.tax.nw.nobad[,189:194])
head(otu.tax.nw.nobad[,1:5])
row.names(otu.tax.nw.nobad) <- otu.tax.nw.nobad[,1] #set first row as rownames
head(otu.tax.nw.nobad[,1:5])
otu.tax.nw.nobad <- otu.tax.nw.nobad[,-1] #remove first row, extraneous Otu column
head(otu.tax.nw.nobad[,1:5])
dim(otu.tax.nw.nobad) #2196 193
head(otu.tax.nw.nobad[,185:193])

#Split again
dim(otu.tax.nw.nobad) #2196 193
head(otu.tax.nw.nobad[,185:193])
otu.notax.nw.nobad <- otu.tax.nw.nobad[,1:191] #take rows 1-193 to make new dataframe otu.notax.nw (192 is delete column, 193 is taxonomy column)
head(otu.notax.nw.nobad[,1:5])
head(otu.notax.nw.nobad[,185:191])
dim(otu.notax.nw.nobad) #2196 191
class(otu.notax.nw.nobad) #dataframe
otu.notax.nw.nobad <- as.matrix(otu.notax.nw.nobad) #turn otu.notax.nw.nobad into a matrix class
class(otu.notax.nw.nobad) #matrix

#Create OTU table phyloseq object
OTUnw.nobad = otu_table(otu.notax.nw.nobad, taxa_are_rows=TRUE)
head(OTUnw.nobad)
dim(OTUnw.nobad) #2196 191
class(OTUnw.nobad)

dim(otu.tax.nw.nobad) #2196 193
head(otu.tax.nw.nobad[,185:193])
tax.levels.nw.nobad <- separate(data = otu.tax.nw.nobad, 
                                col = Taxonomy, 
                                into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#Separate Taxonomy column into 7 separate columns labeled Kingdom to Species
head(tax.levels.nw.nobad) #notice that Species column is blank
dim(tax.levels.nw.nobad) #2196 199
head(tax.levels.nw.nobad[,192:199])
tax.only.nw.nobad <- tax.levels.nw.nobad[,193:198] #keep only taxonomy columns Kingdom to Genus
head(tax.only.nw.nobad)
dim(tax.only.nw.nobad) #2196 6
class(tax.only.nw.nobad) #data.frame
tax.m.nw.nobad <- as.matrix(tax.only.nw.nobad)
class(tax.m.nw.nobad) #matrix
head(tax.m.nw.nobad)

#Create TAX taxonomy table phyloseq object
TAXnw.nobad = tax_table(tax.m.nw.nobad)
head(TAXnw.nobad)
dim(TAXnw.nobad) #2196 6
class(TAXnw.nobad)

#Create phyloseq object containing taxonomy, metadata, and otu table
phyloseq.nw.nobad <- phyloseq(OTUnw.nobad, SAMnw.nobad, TAXnw.nobad)
phyloseq.nw.nobad #view phyloseq object
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 2196 taxa and 191 samples ]
#sample_data() Sample Data:       [ 191 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 2196 taxa by 6 taxonomic ranks ]

save(phyloseq.nw.nobad, file="SRD129.phyloseq.nw.nobad3.RData")

##############################################################################################################################

############# ADONIS TEST FOR NASAL no bad samples, BB and control only, no DNEG 12, DNEG6##############

adonis.nw.nobad <- as(sample_data(phyloseq.nw.nobad), "data.frame")
class(adonis.nw.nobad) #data.frame
dist.nw.nobad <- distance(phyloseq.nw.nobad, method="bray")
set.seed(1)
full.nw.nobad <- adonis(dist.nw.nobad~Day*Treatment, data=adonis.nw.nobad, permutations=9999)
full.nw.nobad

#Call:
#adonis(formula = dist.nw.nobad ~ Day * Treatment, data = adonis.nw.nobad,      permutations = 9999) 

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Day             9    8.9290 0.99211 12.8844 0.36405  1e-04 ***
#Treatment       1    0.6238 0.62377  8.1008 0.02543  1e-04 ***
#Day:Treatment   9    1.8071 0.20079  2.6077 0.07368  1e-04 ***
#Residuals     171   13.1672 0.07700         0.53684           
#Total         190   24.5272                 1.00000 



################PCoA FOR NASAL################

#Distance calculation
#dist.all <- distance(phyloseq.all, "bray") #this didn't work. Got the following error message:
#Error in (function (classes, fdef, mtable)  : 
#unable to find an inherited method for function ‘distance’ for signature ‘"phyloseq", "character"’

#Googled the error message and was given this solution https://github.com/joey711/phyloseq/issues/399:
bray_pcoa_nw <- phyloseq::distance(phyloseq.nw, "bray")

#Perform an ordination on phyloseq data
#PCoA
pcoa.nw <- ordinate(phyloseq.nw, "PCoA", distance=bray_pcoa_nw)

#Reorder Day order
levels(sample_data(phyloseq.nw)$Day) #Output: "D0"     "D1"     "D10"    "D14"    "D21"    "D28"    "D3"     "D36"    "D42"    "D7"     "DNEG12" "DNEG6" 
#looking at SAM aka sample_data Day order
sample_data(phyloseq.nw)$Day <- factor(sample_data(phyloseq.nw)$Day, levels=c("DNEG12", "DNEG6", "D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
#reorder days to appropriate order

#PCoA plot by Day Nasal
pcoa.plot.nw.day <- plot_ordination(phyloseq.nw, pcoa.nw, color="Day") +
  stat_ellipse(type="t") +
  theme(axis.text.x=element_text(color = 'black', size = 13, angle = 45),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Principal Coordinate Analysis plot of nasal samples clustered by day')+  
  theme(plot.title = element_text(size = 17), legend.text = element_text(size=12), legend.title = element_text(size=13), strip.text.x = element_text(size = 10)) +
  scale_color_manual(values = tol12qualitative) #see 12 color palette info below
pcoa.plot.nw.day

#12 color palette from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")

#PCoA by Treatment
pcoa.plot.nw.treatment <- plot_ordination(phyloseq.nw, pcoa.nw, color="Treatment") +
  stat_ellipse(type="t") +
  theme(axis.text.x=element_text(color = 'black', size = 13, angle = 45),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Principal Coordinate Analysis plot of nasal samples clustered by treatment')+  
  theme(plot.title = element_text(size = 17), legend.text = element_text(size=12), legend.title = element_text(size=13), strip.text.x = element_text(size = 10))
pcoa.plot.nw.treatment

#PCoA by treatment and day
pcoa.plot.nw.trtday <- plot_ordination(phyloseq.nw, pcoa.nw, color="All") +
  stat_ellipse(type="t") +
  theme(axis.text.x=element_text(color = 'black', size = 13, angle = 45),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Principal Coordinate Analysis plot of nasal samples clustered by treatment and day')+  
  theme(plot.title = element_text(size = 17), legend.text = element_text(size=12), legend.title = element_text(size=13), strip.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox"))
pcoa.plot.nw.trtday




################PCoA FOR TONSIL################

#Distance calculation
#dist.all <- distance(phyloseq.all, "bray") #this didn't work. Got the following error message:
#Error in (function (classes, fdef, mtable)  : 
#unable to find an inherited method for function ‘distance’ for signature ‘"phyloseq", "character"’

#Googled the error message and was given this solution https://github.com/joey711/phyloseq/issues/399:
bray_pcoa_tt <- phyloseq::distance(phyloseq.tt, "bray")

#Perform an ordination on phyloseq data
#PCoA
pcoa.tt <- ordinate(phyloseq.tt, "PCoA", distance=bray_pcoa_tt)

#Reorder Day order
levels(sample_data(phyloseq.tt)$Day) #Output: "D0"     "D1"     "D10"    "D14"    "D21"    "D28"    "D3"     "D36"    "D42"    "D7"     "DNEG12" "DNEG6" 
#looking at SAM aka sample_data Day order
sample_data(phyloseq.tt)$Day <- factor(sample_data(phyloseq.tt)$Day, levels=c("DNEG12", "DNEG6", "D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
#reorder days to appropriate order

#PCoA plot by Day
pcoa.plot.tt.day <- plot_ordination(phyloseq.tt, pcoa.tt, color="Day") +
  stat_ellipse(type="t") +
  theme(axis.text.x=element_text(color = 'black', size = 13, angle = 45),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Principal Coordinate Analysis plot of tonsil samples clustered by day')+  
  theme(plot.title = element_text(size = 17), legend.text = element_text(size=12), legend.title = element_text(size=13), strip.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox"))
pcoa.plot.tt.day

#PCoA by Treatment
pcoa.plot.tt.treatment <- plot_ordination(phyloseq.tt, pcoa.tt, color="Treatment") +
  stat_ellipse(type="t") +
  theme(axis.text.x=element_text(color = 'black', size = 13, angle = 45),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Principal Coordinate Analysis plot of tonsil samples clustered by treatment')+  
  theme(plot.title = element_text(size = 17), legend.text = element_text(size=12), legend.title = element_text(size=13), strip.text.x = element_text(size = 10))
pcoa.plot.tt.treatment

#PCoA by treatment and day
pcoa.plot.tt.trtday <- plot_ordination(phyloseq.tt, pcoa.tt, color="All") +
  stat_ellipse(type="t") +
  theme(axis.text.x=element_text(color = 'black', size = 13, angle = 45),
        axis.text.y=element_text(color = 'black', size=13), 
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15))+ ggtitle('Principal Coordinate Analysis plot of tonsil samples clustered by treatment and day')+  
  theme(plot.title = element_text(size = 17), legend.text = element_text(size=12), legend.title = element_text(size=13), strip.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox"))
pcoa.plot.tt.trtday


######################## Jules Zone ########################
library(vegan)

nw2.meta <- data.frame(phyloseq.nw2@sam_data) #make phyloseq.nw2 sam_data into dataframe
nw2.otu <- data.frame(t(phyloseq.nw2@otu_table)) #make phyloseq.nw2 otu_table into dataframe
class(nw2.meta) #data.frame
rownames(nw2.meta) == row.names(nw2.otu) #make sure rownames between nw2.meta and nw2.otu match exactly. Yes they do.

#Beta-diversity stats with pairwise adonis
nw.adon <- pairwise.adonis(nw2.otu, nw2.meta$All) #run pairwise adonis on otu table and All column of meta data
#nw.adon contains all the pairwise comparisons
nw.adon$pairs #list all comparisons in the "pairs" column
goodcomps <- c(grep('D0 [A-Za-z]+ vs D0 [A-Za-z]+', nw.adon$pairs),
  grep('D4 [A-Za-z]+ vs D4 [A-Za-z]+', nw.adon$pairs),
  grep('D7 [A-Za-z]+ vs D7 [A-Za-z]+', nw.adon$pairs),
  grep('D11 [A-Za-z]+ vs D11 [A-Za-z]+', nw.adon$pairs),
  grep('D14 [A-Za-z]+ vs D14 [A-Za-z]+', nw.adon$pairs))
#use regular expressions
#[A-Za-z] match all capital, lowercase letters
#+ match a whole word and not just one letter (if you didn't have "+")
#c creates the vector, lump all pairs of specific interest groups together
#make vector to combine all the pairwise comparisons you're interested in (same day, different treatment group)

nw.adon.good <- nw.adon[goodcomps,] #rename this vector

#Alpha-diversity stats with pairwise wilcox test
nw2.meta$shannon <- diversity(nw2.otu) #diversity is vegan function, default index is shannon
#added a shannon index column in nw2.meta


meta0 <- nw2.meta[nw2.meta$Day == 'D0',] #take day 0 rows and create meta0 object
pairwise.wilcox.test(meta0$shannon, meta0$All) #run pairwise wilcox test on meta0 (day 0 comparisons)

typeof(nw2.meta) #list... need to make into table
dim(nw2.meta) #207 6
head(nw2.meta[,1:6]) #shannon column at the end
pairwise.wilcox.test(nw2.meta$shannon, nw2.meta$All, p.adjust.method = 'none')

#Alpha-diversity box whisker plot (based on "All" column)
#x-axis day
#legend: treatments
#make for tonsil, nasal separate box whisker plots

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
