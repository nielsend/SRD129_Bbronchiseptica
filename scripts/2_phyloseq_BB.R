#########################################
#SRD129 16S - Creating phyloseq objects

#Purpose: Create phyloseq objects to be used to calculate alpha and beta diversity measures for nasal samples in BB and Control groups. 
#This section will also use the adonis function to determine the effect of time and treatment on the community structure of nasal microbiota.

#Files needed:
#OTU table: SRD129BBabundsingleton2000OTUtable.csv
#Metadata: SRD129BBmetadata.csv

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)

#Read files
otu <- read.csv("./SRD129BBabundsingleton7395OTUtable.csv", row.names=1) #Set column 1 as row names
meta <-read.csv("./metadata.csv")
colnames(meta) [1] <- "Sample" #Rename first column of 'meta' to "Sample"

dim(otu) #Check dimensions of 'otu' = 1636 rows and 178 columns
head(otu[,1:5]) #Check the first part of 'otu' table
head(otu[,170:178]) #Check last part of 'otu' table
dim(meta) #Determine the dimensions of 'meta' dataframe. It has 190 rows and 5 columns
head(meta[,1:5]) #Check first part of 'meta' table

#Remove taxonomy from 'otu' DANIEL YOU ARE HERE
tax <- otu[,(177:178)] #Copy column 177 taxonomy and 178 to keep the row names from 'otu' to 'tax'
head(tax)
colnames(tax)[1] <- "delete" #Rename column 1 of 'tax' (formerly column 177) as "delete" which will be deleted later
head(tax)

#Modify 'otu' with only OTU count data
head(otu)
dim(otu)
otu <- otu[,-178] #Remove column 178 taxonomy in 'otu' to have only OTU data
dim(otu) #Dimensions of 'otu' show 1636 rows 177 columns
head(otu[,170:177])

#Transpose 'otu' to match format of 'meta'
otu.trans <- t(otu)
dim(otu.trans) #177 1636
#Now rownames in 'otu.trans' are sample names, columns are OTUs
head(otu.trans[,1:5])
head(meta) #Row names are numbered, but we want sample names as row names
rownames(meta) <- meta$Sample #Set names in "Sample" as rownames in 'meta'
colnames(meta)[1] <- "Sample"

head(meta)
meta <- meta[,-1]
dim(meta) #190   4
class(meta) #The type of class that 'meta' is is dataframe
class(otu) #data.frame

#Merge 'otu' and 'meta' data frames
otu.meta <- merge(meta, otu.trans, by.x=0, by.y=0) #Merge by names of the columns that are common to both x and y (columns with common names between the two data sets)
#by.x=0 means match by rownames in 'meta'; y=0 means match by rownames in 'otu.trans'
dim(otu.meta)
head(otu.meta[,1:10])
class(otu.meta) #Check class type of 'otu.meta'. It should be a dataframe.

#Added an "All" column (combines 'Treatment' and 'Day' values) in new 'otu.meta2' dataframe
otu.meta2<- cbind(otu.meta) #Make second copy of 'otu.meta' to use to include an "All" column
otu.meta2$All <- with(otu.meta2, paste0(Day, sep="_", Treatment)) #Combine "Day" and "Treatment" columns into an "All" column
dim(otu.meta2) #176 1642
head(otu.meta2[,1636:1642]) #Check the first part of the end of 'otu.meta2' dataframe
rownames(otu.meta2) <- otu.meta2$Row.names #Set "Row.names" as rownames
otu.meta2 <- otu.meta2[,-1]
head(otu.meta2[,1:10])
dim(otu.meta2) #176 1641
head(otu.meta2[,1631:1641])
otu.meta2<- otu.meta2[,c(1:4,1641,5:1640)] #Reorder columns to have "All" column after "Treatment" column
write.csv(otu.meta2, file="SRD129BBabundsingleton7395.otu.meta.csv")

##DANIEL YOU ARE HERE


#Pull out metadata from 'otu.meta2' dataframe
head(otu.meta2[,1:10])
dim(otu.meta2) #176 1641
otu.meta3 <- otu.meta2[,1:5] #Take columns 1-5 to make 'otu.meta3'
dim(otu.meta3) #176 5
head(otu.meta3)

#Create SAM metadata table phyloseq object
SAM = sample_data(otu.meta3, errorIfNULL = TRUE)
head(SAM)
dim(SAM) #176 5

#Pull out OTU data from 'otu.meta2' dataframe
head(otu.meta2[,1631:1641])
otu.meta4 <- otu.meta2[,c(6:1641)] #Select "OTU" columns to create 'otu.meta4' dataframe
head(otu.meta4[,1:10])
dim(otu.meta4) #176 1636
otu.meta4.trans <- t(otu.meta4) #Transpose 'otu.meta4' to have OTUs as rownames, sample names as column names
head(otu.meta4.trans[,1:10])
dim(otu.meta4.trans) #1636 176

#Merge 'tax' back into 'otu.meta4.trans' for correct format and taxons
head(tax)
otu.tax <- merge(otu.meta4.trans, tax, by=0) #Merge by rownames aka OTU rownames
dim(otu.tax) #1636 179
head(otu.tax[,171:179])
head(otu.tax[,1:5])
row.names(otu.tax) <- otu.tax[,1] #Set first row as rownames
head(otu.tax[,1:5])
otu.tax <- otu.tax[,-1] #Remove first row, extraneous OTU column
head(otu.tax[,1:5])
dim(otu.tax) #1636 178
head(otu.tax[171:178])

#Split 'otu.tax.nw2' again
otu.notax <- otu.tax[,1:176] #Take rows 1-176 to make new dataframe 'otu.notax' (177 is delete column, 178 is taxonomy column)
head(otu.notax[,1:5])
dim(otu.notax) #1636 176
class(otu.notax) #dataframe
otu.notax <- as.matrix(otu.notax) #Turn 'otu.notax' into a matrix class
class(otu.notax) #matrix, array

#Create OTU table phyloseq object
OTU = otu_table(otu.notax, taxa_are_rows = TRUE)
dim(OTU) #1636 176
class(OTU)
dim(otu.tax) #1636 178
tax.levels <- separate(data = otu.tax,
                       col = Taxonomy,
                       into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")

#Separate Taxonomy column into 7 separate columns labeled "Kingdom", "Phylum", "Class", etc.
head(tax.levels) #Notice that "Species" column is blank
dim(tax.levels) #1636 184
head(tax.levels[,176:184])
tax.only <- tax.levels[,178:183] #Keep only taxonomy columns from "Kingdom" up to "Genus"
head(tax.only)
dim(tax.only) #1636 6
class(tax.only) #data.frame
tax.m <- as.matrix(tax.only)
class(tax.m) #matrix
head(tax.m)

#Create TAX taxonomy table phyloseq object
TAX = tax_table(tax.m)
head(TAX)
dim(TAX) #1636 6
head(TAX)
class(TAX)

#Create phyloseq object 'phyloseqFlu' containing taxonomy, metadata, and OTU table
phyloseqbb <- phyloseq(OTU, SAM, TAX)
phyloseqbb #view phyloseq object
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1636 taxa and 176 samples ]
#sample_data() Sample Data:       [ 176 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 1636 taxa by 6 taxonomic ranks ]

save(phyloseqbb, file="SRD129bb.phyloseq.RData")

#Run adonis function to determine effect of time and treatment on structure of nasal microbiota
adonis.bb <- as(sample_data(phyloseqbb), "data.frame")
class(adonis.bb) #data.frame

#distance function to run distance calculations
dist.bb <- distance(phyloseqbb, method="bray") #Distance calculation using Bray-Curtis
set.seed(1) #Use set.seed function when running simulations to ensure all results are reproducible
full.bb <- adonis(dist.bb~Day*Treatment, data=adonis.bb, permutations=9999)
full.bb #Display results





#switched Treatment and Day order
full.bb_2 <- adonis(dist.bb~Treatment*Day, data=adonis.bb, permutations=9999)
full.bb_2 #Display results



#If distance function is giving an error message below:
#"Error: x should be a data.frame, data.table, tbl, tbl_df, array, or matrix."
#Use vegdist function from vegan package to run distance calculations instead of the distance function
#and use those calculations to run through adonis test.
#vegdist requires that phyloseq object's OTU table has OTUs listed in the columns and sample names listed in rows.
#Also, remove any OTUs with taxa_sums = 0 or non-numeric values. For example, this command can help remove OTUs with taxa_sums = 0:
#OTU <- prune_taxa(taxa_sums(<yourOTUtable>) > 0, <yourOTUtable>)
#If you create a separate phyloseq object with this specific OTU table setup,
#you should be able to run the vegdist function without any errors and use the output to run through adonis function.

head(otu.meta4) #Sample names are listed in rows and OTUs are listed in columns in 'otu.meta4'
OTU.2 = otu_table(otu.meta4, taxa_are_rows = TRUE)
OTU.2.distance <- vegdist(OTU.2, method = "bray")
OTU.2.distance.full <- adonis(OTU.2.distance~Day*Treatment, data = adonis.bb, permutations=9999)
OTU.2.distance.full
#Call:
#adonis(formula = OTU.2.distance ~ Day * Treatment, data = adonis.bb,      permutations = 9999) 

# $aov.tab
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
#   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Day             8    8.0310 1.00387 13.6468 0.36619  1e-04 ***
#   Treatment       1    0.6411 0.64108  8.7150 0.02923  1e-04 ***
#   Day:Treatment   8    1.6365 0.20456  2.7809 0.07462  1e-04 ***
#   Residuals     158   11.6226 0.07356         0.52996           
# Total         175   21.9312                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $call
# adonis(formula = OTU.2.distance ~ Day * Treatment, data = adonis.bb, 
#        permutations = 9999)
# 
# $coefficients
# NULL

head(otu.meta4) #Sample names are listed in rows and OTUs are listed in columns in 'otu.meta4'
OTU.2 = otu_table(otu.meta4, taxa_are_rows = TRUE)
OTU.2.distance <- vegdist(OTU.2, method = "bray")
OTU.2.distance.full <- adonis(OTU.2.distance~Treatment*Day, data = adonis.bb, permutations=9999)
OTU.2.distance.full

#adonis(OTU.2.distance~Treatment*Day, data = adonis.bb, permutations=9999)
# $aov.tab
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Treatment       1    0.6230 0.62297  8.4687 0.02841  1e-04 ***
#   Day             8    8.0491 1.00614 13.6776 0.36702  1e-04 ***
#   Treatment:Day   8    1.6365 0.20456  2.7809 0.07462  1e-04 ***
#   Residuals     158   11.6226 0.07356         0.52996           
# Total         175   21.9312                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $call
# adonis(formula = OTU.2.distance ~ Treatment * Day, data = adonis.bb, 
#        permutations = 9999)
# 
# $coefficients
# NULL