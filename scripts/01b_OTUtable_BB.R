#########################################
#SRD129 16S - OTU table
#By Mou, KT

#Purpose: Create OTU table with R using specific files generated from mothur

#Files needed:
#SRD129BB.outsingletons.abund.opti_mcc.0.03.cons.taxonomy
#SRD129BB.outsingletons.abund.taxonomy.csv
#SRD129BB.outsingletons.abund.opti_mcc.0.03.subsample.shared
#SRD129BB.outsingletons.abund.subsample.shared.csv

#The output files from mothur needed in R for this section and subsequent sections, aside from fasta files, are text files that can be saved as csv for ease of use in R.

#Load library package
library(tidyverse)

#To start creating OTU table, edit taxonomy file
#Made a copy of "SRD129BB.outsingletons.abund.opti_mcc.0.03.cons.taxonomy" file, opened the file in a spreadsheet editor and removed the "Size" column. Renamed the copy to "SRD129BB.outsingletons.abund.taxonomy.csv".
taxonomy <- read.csv("./data/SRD129BB.outsingletons.abund.taxonomy.csv")  #Import this csv file from working directory using "read.csv" function
taxonomy$Taxonomy <- gsub('*\\(.*?\\) *', '', taxonomy$Taxonomy) #Substitute (100) and variations of that with nothing ('') from the Taxonomy column
taxonomy[1:6,2] #Show column number 2, rows 1 through 6 in 'taxonomy' dataframe
write.csv(taxonomy, file = "SRD129BBabundsingleton2000taxonomy.csv")

#Edit subsample.shared file
#Save "SRD129BB.outsingletons.abund.opti_mcc.0.03.subsample.shared" file generated from mothur as a csv file in a spreadsheet editor. Named file as "SRD129.outsingletons.abund.subsample.shared.csv".
shared <- read.csv("./data/SRD129BB.outsingletons.abund.subsample.shared.csv")
head(shared) #Check on the first part of 'shared' dataframe
shared[1:6,1]
shared <- t(shared) #Transpose 'shared'
head(shared)
shared[1:6,1]
write.csv(shared, file = 'SRD129BBabundsingleton2000shared.csv') #Open this csv file in a spreadsheet editor and remove the "V*", "label", and "numOtus" rows
#numOtus = 1308

#Combine shared and taxonomy files to create OTU table
shared <- read.csv("./data/SRD129BBabundsingleton2000shared.csv")
colnames(shared) [1] <- "OTU" #Rename first column of 'shared' to "OTU"
taxonomy <- read.csv("./data/SRD129BBabundsingleton2000taxonomy.csv")
OTUtable <- merge(shared, taxonomy, by.x ="OTU", by.y = "OTU") #Merge 'shared' and 'taxonomy' objects by OTU
head(OTUtable)
nrow(OTUtable) #Count number of rows in 'OTUtable'
ncol(OTUtable) #Count number of columns in 'OTUtable'
write.csv(OTUtable, file= "SRD129BBabundsingleton2000OTUtable.csv") 
#Open this csv file in a spreadsheet editor and remove the first column (numbered rows) and column X (second to last column)

#Check OTU table
OTUtable <-read.csv("./data/SRD129BBabundsingleton2000OTUtable.csv", stringsAsFactors = FALSE)
head(OTUtable) #Check the first part of 'OTUtable' to make sure it looks ok
