#########################################
#SRD129 16S - OTU table


#Purpose: Create OTU table with R using specific files generated from mothur

#Files needed:
#SRD129BB.outsingletons.abund.opti_mcc.0.03.cons.taxonomy
#SRD129BB.outsingletons.abund.taxonomy.csv
#SRD129BB.outsingletons.abund.opti_mcc.0.03.subsample.shared
#SRD129BB.outsingletons.abund.subsample.shared.csv

#The output files from mothur that are needed in R for this section and subsequent sections, 
#aside from fasta files, are text files that can be saved as csv for ease of use in R.

#Load library package
library(tidyverse)

#To start creating OTU table, edit taxonomy file
#Use "SRD129BB.outsingletons.abund.opti_mcc.0.03.cons.taxonomy" file
taxonomy <- read_tsv("SRD129BB.outsingletons.abund.opti_mcc.0.03.cons.taxonomy")  #Import this csv file from working directory using "read.csv" function
taxonomy$Taxonomy <- gsub('*\\(.*?\\) *', '', taxonomy$Taxonomy) #Substitute (100) and variations of that with nothing ('') from the Taxonomy column
taxonomy[1:6,2] #Show column number 2, rows 1 through 6 in 'taxonomy' dataframe
write.csv(taxonomy, file = "SRD129BBabundsingleton7395taxonomy.csv")

#Edit subsample.shared file
#Use "SRD129BB.outsingletons.abund.opti_mcc.0.03.subsample.shared" file generated from mothur 
shared <- read_tsv("./SRD129BB.outsingletons.abund.opti_mcc.0.03.subsample.shared")
head(shared) #Check on the first part of 'shared' dataframe
shared[1:6,1]
shared <- t(shared) #Transpose 'shared'
head(shared)
shared[1:6,1]
colnames(shared) <- shared[2,] # set colnames as the group names
shared <- shared[-(1:3),] # remove the "V*", "label", and "numOtus" rows
write.csv(shared, file = 'SRD129BBabundsingleton7395shared.csv') #save DF
dim(shared)
#numOtus = 1636 dwn

#Combine shared and taxonomy files to create OTU table
shared <- read.csv("./SRD129BBabundsingleton7395shared.csv")
colnames(shared) [1] <- "OTU" #Rename first column of 'shared' to "OTU"
taxonomy <- read.csv("./SRD129BBabundsingleton7395taxonomy.csv")
OTUtable <- merge(shared, taxonomy, by.x ="OTU", by.y = "OTU") #Merge 'shared' and 'taxonomy' objects by OTU
head(OTUtable)
nrow(OTUtable) #Count number of rows in 'OTUtable'
ncol(OTUtable) #Count number of columns in 'OTUtable'
#remove the first column (numbered rows) and column X (second to last column)
OTUtable <- OTUtable %>% select(-`X`) #%>%  select(-`Group`) %>% 
write.csv(OTUtable, file = "SRD129BBabundsingleton7395OTUtable.csv", row.names=FALSE)
OTUtable <-read.csv("SRD129BBabundsingleton7395OTUtable.csv", stringsAsFactors = FALSE)
head(OTUtable) #Check the first part of 'OTUtable' to make sure it looks ok


View(OTUtable)

OTUtable$SRD129_P141_N_D3
NotInCounts <- OTUtable %>%  select(c("OTU", "SRD129_P141_N_D3", "SRD129_P133_N_D1"))
NotInCounts %>% filter(OTU=="Otu0032")


Bquestion <- OTUtable %>%  filter(OTU=="Otu0032")
Bquestion <- t(Bquestion)
Bquestion <- Bquestion[-1,]
Bquestion <-as.data.frame(Bquestion)
Bquestion$Sample <- rownames(Bquestion)
write.csv(Bquestion, "./Bquestion.csv")
colnames(Bquestion) <- c("OTU", "Sample")
ggplot(Bquestion, aes(x=Sample, y=as.numeric(OTU))) + geom_bar(stat="identity") +  theme(axis.text.x=element_text(size=12, angle=45, hjust=1))

