#######################################################################
#SRD129 Bordetella bronchiseptica KM22 nasal count data
#Kathy Mou

#NOTES: 
#This code graphs the KM22 CFU counts from nasal wash samples, as box and whisker plots 

#Clear workspace and load necessary packages
rm(list=ls())

#Set working directory on Mac desktop:
setwd("~/Desktop/Microbiome/Projects/SRD129/SRD129_2000singletons")
setwd("~/Desktop/SRD129/SRD129_2000singletons")

#Load library packages
library(ggplot2)


#Load image file
load("SRD129_BBKM22Count.RData")



#Save image file
save.image(file="SRD129_BBKM22Count.RData")

#######################################################################

bb<- read.csv('BBKM22CFUCount.csv')
bb$Day = factor(bb$Day, levels=c("DNEG12", "DNEG6", "D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
factor(bb$Day) #Levels: DNEG12 DNEG6 D0 D1 D3 D7 D10 D14 D21 D28 D36 D42
bb.1<- bb[!grepl("DNEG12", bb$Day),]
bb.1<- bb.1[!grepl("DNEG6", bb.1$Day),]
(b <- ggplot(bb.1, aes(Day, CFUpermL)))
b2 <- b+ geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle=45, hjust=1)) + ylab("CFU per mL") +
      scale_y_continuous(trans=log10_trans())
b2
