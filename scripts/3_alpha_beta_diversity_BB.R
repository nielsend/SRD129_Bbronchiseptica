#########################################
#SRD129 16S - Alpha and Beta diversity
#By Mou, KT

#Purpose: This code generates non-metric multidimensional scaling ordination based on Bray-Curtis dissimilarities to create NMDS plots, and runs pairwise.adonis function to identify any significant differences in bacterial composition between treatment groups on a given day.

#Files needed:
#phyloseqbb

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)
library(scales)
library(RColorBrewer)
library(philentropy)
library(cowplot)
library('wesanderson')

#Load this function from the funfuns R package (https://github.com/Jtrachsel/funfuns)
NMDS_ellipse <- function(metadata, OTU_table, grouping_set,
                         distance_method = 'bray',
                         rand_seed = 77777,
                         MDS_trymax = 1000,
                         autotransform = FALSE,
                         wascores = TRUE,
                         expand = FALSE){
  require(vegan)
  require(tidyr)
  
  if (grouping_set %in% colnames(metadata)){
    if (all(rownames(metadata) == rownames(OTU_table))){
      
      set.seed(rand_seed)
      generic_MDS <- metaMDS(OTU_table, k = 2,
                             trymax = MDS_trymax,
                             autotransform = autotransform,
                             distance = distance_method,
                             wascores = wascores,
                             expand = expand)
      
      stress <- generic_MDS$stress
      nmds_points <- as.data.frame(generic_MDS$points)
      metadata <- metadata[match(rownames(generic_MDS$points), rownames(metadata)),]
      metadata <- as.data.frame(metadata) # weird things were happening when a grouped tibble was used as metadata...
      metanmds <- cbind(metadata, nmds_points)
      # browser()
      nmds.mean <- aggregate(metanmds[,grep("MDS", colnames(metanmds))], list(group=metanmds[[grouping_set]]), mean)
      metanmds[[grouping_set]] <- factor(metanmds[[grouping_set]]) # this 'set' needs to be passed in from function
      
      #check to make sure at least 3 obs for each grouping_set
      
      numobs <- metanmds %>% group_by(!!grouping_set) %>% summarise(n=n())
      if (min(numobs$n) >= 3){
        ord <- ordiellipse(generic_MDS, metanmds[[grouping_set]], label = TRUE, conf = .95, kind = 'se', draw = 'none')
        
        df_ell <- data.frame()
        for (d in levels(metanmds[[grouping_set]])){
          df_ell <- rbind(df_ell, cbind(as.data.frame(with(metanmds[metanmds[[grouping_set]] == d,],
                                                           vegan:::veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
        }
        
        
        
        # this loop assigns the group centroid X coordinates to each sample, there is probably a better way...
        
        metanmds$centroidX <- NA
        metanmds$centroidY <- NA
        
        
        for (level in levels(metanmds[[grouping_set]])){
          metanmds[metanmds[[grouping_set]] == level,]$centroidX <- nmds.mean$MDS1[nmds.mean$group == level]
          metanmds[metanmds[[grouping_set]] == level,]$centroidY <- nmds.mean$MDS2[nmds.mean$group == level]
          
          
        }
        print(paste('Ordination stress:', stress, sep = ' '))
        return(list(metanmds, df_ell, generic_MDS))
        
      } else {
        warning('One of your groups in "grouping_set" has less than 3 observations, cannot generate elipses')
        df_ell <- data.frame()
        return(list(metanmds, df_ell, generic_MDS))}
      
    } else {
      stop('The rownames for your OTU table and metadata do not match.')
    }
    
  }else {
    stop('The grouping set column you have provided in not in your metadata.')
  }
  
  
  
}


#Load the following pairwise.adonis function, taken from https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
}

#######################################################################

#Setting up 'phyloseqbb' into dataframes for beta diversity NMDS calculation
bb.sam <- data.frame(phyloseqbb@sam_data) #Make 'phyloseqbb sam_data' into dataframe. Carry over phyloseqbb object from 02_phyloseq_BB.R
bb.otu <- data.frame(t(phyloseqbb@otu_table)) #Make 'phyloseqbb otu_table' into dataframe
class(bb.sam) #data.frame
rownames(bb.sam) == row.names(bb.otu) #For rows with sums greater than 1 in 'bb.otu', move rows and their respective sum values into "numOTUs" column in 'bb.sam'
bb.sam$numOTUS <- rowSums(bb.otu > 1)
head(bb.sam)

#NMDS calculation
bb.otu[1:10,1:10]
bb.NMDS <- NMDS_ellipse(bb.sam, bb.otu, grouping_set = 'All')
#Output:
#Result: "Ordination stress: 0.16936675872559"

#Separate meta data and ellipse data to two lists to make NMDS plot
bb.NMDS.2 <- bb.NMDS[[1]]
#'bb.NMDS.2' has meta data + MDS calculations. Select this 1st list of 'bb.NMDS' using double brackets
bb.df_ell.2 <- bb.NMDS[[2]]
#'bb.df_ell.2' is accessing 2nd list from 'bb.NMDS' that has ellipse calculations
#Need two separate lists for making NMDS plot
bb.df_ell.2$group
head(bb.df_ell.2)

#Create "Day" and "Treatment" columns within 'bb.df_ell.2' for faceting purposes
bb.df_ell.2$Day <- sub('_[A-Za-z]+', '\\1', bb.df_ell.2$group) #Created "Day" column, '\\1' returns the first part of the regular expression D[A-z0-9]+ from 'bb.df_ell.2$group'
bb.df_ell.2$Treatment <- sub('D[A-z0-9]+_', '\\2', bb.df_ell.2$group) #Create "Treatment" column, '\\2' returns the second part of the sub expression ([A-Za-z]+) from 'bb.df_ell.2$group'
head(bb.df_ell.2)

#Restructure "Day" level order for 'bb.df_ell.2' and 'bb.NMDS.2'
bb.NMDS.2$Day = factor(bb.NMDS.2$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D36", "D42"))
bb.df_ell.2$Day = factor(bb.df_ell.2$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D36", "D42"))
levels(bb.df_ell.2$Day) # [1] "D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D36" "D42"
levels(bb.NMDS.2$Day) # [1] "D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D36" "D42"
View(bb.NMDS.2)
#Creating plot from NMDS calculations
bbNMDSplot <- ggplot(data=bb.NMDS.2, aes(x=MDS1, y=MDS2, color=Treatment)) + geom_point() +
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = 0.5) +
  geom_path(data=bb.df_ell.2, aes(x=NMDS1, y=NMDS2, color=Treatment, group=group)) +
  facet_wrap(~Day, scales = 'free') +
  #scale_color_brewer(palette="Dark2") +
  theme_gray(base_size = 10) +
  theme(strip.text.x = element_text(size=15)) +
  labs(color="Treatment group") + labs(caption = 'Ordination stress = 0.169') +
  theme(legend.text = element_text(size=13),
        legend.title = element_text(size=13))
bbNMDSplot

#Save 'bbNMDSplot' as a .tiff for publication, 500dpi
ggsave("bbNMDSplot1.tiff", plot=bbNMDSplot, width = 10, height = 6, dpi = 500, units =c("in"))

#Using pairwise.adonis function
bb.adon <- pairwise.adonis(bb.otu, bb.sam$All, sim.method = 'bray', p.adjust.m = 'bonferroni')
#Run pairwise.adonis on 'bb.otu' OTU table and "All" column of 'bb.sam' dataframe
bb.adon$pairs #List all comparisons in the "pairs" column of 'bb.adon'
goodcomps.bb.adon <- c(grep('D0_[A-Za-z]+ vs D0_[A-Za-z]+', bb.adon$pairs),
                       grep('D1_[A-Za-z]+ vs D1_[A-Za-z]+', bb.adon$pairs),
                       grep('D3_[A-Za-z]+ vs D3_[A-Za-z]+', bb.adon$pairs),
                       grep('D7_[A-Za-z]+ vs D7_[A-Za-z]+', bb.adon$pairs),
                       grep('D10_[A-Za-z]+ vs D10_[A-Za-z]+', bb.adon$pairs),
                       grep('D14_[A-Za-z]+ vs D14_[A-Za-z]+', bb.adon$pairs),
                       grep('D21_[A-Za-z]+ vs D21_[A-Za-z]+', bb.adon$pairs),
                       grep('D36_[A-Za-z]+ vs D36_[A-Za-z]+', bb.adon$pairs),
                       grep('D42_[A-Za-z]+ vs D42_[A-Za-z]+', bb.adon$pairs))
# "[A-Za-z]" matches all capital and lowercase letters
# "+" matches a whole word and not just one letter (if you didn't have "+", then it would match by one letter)
# "c" creates the vector, lumps all pairs of specific groups of interest together
# You want to make a vector to combine all the pairwise comparisons you're interested in (same day, different treatment group)
goodcomps.bb.adon
goodcomps.bb.adon.2 <- bb.adon[goodcomps.bb.adon,] #Rename 'goodcomps.bb.adon' vector to 'goodcomps.bb.adon.2'
goodcomps.bb.adon.2
goodcomps.bb.adon.2$p.adjusted <- p.adjust(goodcomps.bb.adon.2$p.value, method = 'fdr') #"p.adjust" function returns a set of p-values adjusted with "fdr" method
goodcomps.bb.adon.2$p.adjusted2 <- round(goodcomps.bb.adon.2$p.adjusted, 4) #Round p-values to 3 decimal points and list in new "p.adjusted2" column
write.csv(goodcomps.bb.adon.2, file='permanova1.csv', row.names=TRUE)
goodcomps.bb.adon.2$p.adjusted2[goodcomps.bb.adon.2$p.adjusted2 > 0.05] <- NA #For all p-values greater than 0.05, replace with "NA"
write.csv(goodcomps.bb.adon.2, file='bb.pairwisecomparisons.csv', row.names=TRUE)

#Alpha diversity: average shannon, invsimpson, numOTUs per "All" subtype
bb.sam$shannon <- diversity(bb.otu) #diversity is vegan function, default index is shannon
#added a shannon index column in nw.meta.nobad
bb.sam$invsimpson <- diversity(bb.otu,index = 'invsimpson')
#use invsimpson since it is easier to understand than simpson
#don't need to "inverse" simpson values
#simpson: low number = higher diversity

levels(sample_data(bb.sam)$Day) #"D0"  "D1"  "D10" "D14" "D21" "D3"  "D36" "D42" "D7"
bb.sam$Day = factor(bb.sam$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D36", "D42"))
levels(sample_data(bb.sam)$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"   
head(bb.sam)

#use underscores instead of spaces
#bb.sam$All <- gsub(" ", "_", bb.sam$All)

#average shannon, invsimpson, numOTUs per "All" subtype
bb.sam.average.shannon.invsimpson.numOTUs <- aggregate(bb.sam[, 6:8], list(bb.sam$All), mean)
print(bb.sam.average.shannon.invsimpson.numOTUs)

# Group.1  numOTUS  shannon invsimpson
# Group.1   numOTUS  shannon invsimpson
# 1        D0_BB  94.60000 2.081025   4.456822
# 2   D0_Control 145.90000 2.597476   6.312613
# 3        D1_BB  97.44444 2.312722   6.005666
# 4   D1_Control 106.10000 2.276661   5.115462
# 5       D10_BB 146.50000 2.707806   7.482703
# 6  D10_Control 132.70000 2.489575   5.186045
# 7       D14_BB 151.50000 2.517312   4.075226
# 8  D14_Control 167.60000 2.888996   6.956794
# 9       D21_BB 187.10000 3.142255   7.980939
# 10 D21_Control 159.50000 2.944002   7.078478
# 11       D3_BB 165.55556 2.977226   7.583899
# 12  D3_Control 160.20000 2.945695   7.385418
# 13      D36_BB 200.70000 3.463784  13.791708
# 14 D36_Control 192.10000 3.208963  13.187077
# 15      D42_BB 260.80000 3.996025  29.422445
# 16 D42_Control 218.30000 3.561202  17.418921
# 17       D7_BB 135.37500 2.420744   5.385052
# 18  D7_Control 116.50000 2.421066   5.456063

bb.pairwise.wilcox.shannon.test <- pairwise.wilcox.test(bb.sam$shannon, bb.sam$All, p.adjust.method = 'none') #Calculate pairwise comparisons by "All" column of the shannon indices in "Shannon" column
print(bb.pairwise.wilcox.shannon.test) #Look at the results of 'bb.pairwise.wilcox.shannon.test'
bb.pairwise.wilcox.invsimpson.test <- pairwise.wilcox.test(bb.sam$invsimpson, bb.sam$All, p.adjust.method = 'none')
print(bb.pairwise.wilcox.invsimpson.test)

#box and whisker plot: shannon
bb.shan.bw <- ggplot(data = bb.sam, aes(x=Treatment, y=shannon, group=All, fill=Treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~Day, scales = 'free_x') +
  scale_y_continuous(name = "Shannon diversity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text.x = element_text(size=10),
        axis.title.y = element_text(size=15))
#free allows each plot to customize scale to the specific data set (no forced scale applied to all plots)
#position = position_dodge2(preserve = 'total') --> fixes the ggplot box width;
#prevents narrow boxes from forming in plot (they become wider)
bb.shan.bw
ggsave("Figure_Shannon.tiff", plot=bb.shan.bw, width = 10, height = 5, dpi = 500, units =c("in"))



#box and whisker plot: inverse simpson
bb.invsimp.bw <- ggplot(data = bb.sam, aes(x=Treatment, y=invsimpson, group=All, fill=Treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~Day) +
  scale_y_continuous(name = "Inverse Simpson diversity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text.x = element_text(size=10),
        axis.title.y = element_text(size=15))
bb.invsimp.bw

#Save 'bb.invsimp.bw' as a .tiff for publication, 500dpi
ggsave("Figure_InvSimpson.tiff", plot=bb.invsimp.bw, width = 10, height = 5, dpi = 500, units =c("in"))

#Combining plots
bb.combalpha <- plot_grid(bb.shan.bw + theme(legend.position = "none"), bb.invsimp.bw, labels = "AUTO", rel_widths = c(45.5,54.5), scale=0.9)
bb.combalpha

#Save 'flu.combalpha' as a .tiff for publication, 500dpi
ggsave("Figure_1.tiff", plot=bb.combalpha, width = 14, height = 8, dpi = 500, units =c("in"))





#sum Total of OTU in each column
sum <- 0
my_matrix <- otu.meta4.trans
for(i in 1:ncol(my_matrix)){
  sum[i] <- sum(my_matrix[,i])
}

#determine unique number of OTU in each column
uniqueOTUinSubsample <- as.data.frame(apply(my_matrix,2,function(x) sum(x > 0)))
uniqueOTUinSubsample <- as.data.frame(uniqueOTUinSubsample)
uniqueOTUinSubsample$group <- rownames(uniqueOTUinSubsample)
colnames(uniqueOTUinSubsample) <- c("numOTU", "group")

#DANIEL YOU ARE HERE
meta2 <- meta
meta2$group <- rownames(meta2)
uniqueOTUinSubsampleTot <- merge(meta2, uniqueOTUinSubsample, by="group")
uniqueOTUinSubsampleTot$Set <- paste(uniqueOTUinSubsampleTot$Day, uniqueOTUinSubsampleTot$Treatment)
View(uniqueOTUinSubsampleTot)

# arrange for graphpad
GPuniqueOTU <- uniqueOTUinSubsampleTot %>% group_by(`Set`)
GPuniqueOTU <- GPuniqueOTU %>% subset(select= - Day)
GPuniqueOTU <- GPuniqueOTU %>% subset(select= - Treatment)
write.csv(GPuniqueOTU, "./GPunique.OTU.csv")


# generate sum per group & day
uniOTUSubSampleMean <- data.frame(aggregate(uniqueOTUinSubsampleTot$numOTU, by=list(Category=uniqueOTUinSubsampleTot$Set), FUN=mean))
uniOTUSubSampleSD <- data.frame(aggregate(uniqueOTUinSubsampleTot$numOTU, by=list(Category=uniqueOTUinSubsampleTot$Set), FUN=sd))
uniOTUSubSampleMean <- merge(uniOTUSubSampleMean, uniOTUSubSampleSD, by="Category")
uniOTUSubSampleMean <- uniOTUSubSampleMean %>% mutate_if(is.numeric, round, .1) #round mean
colnames(uniOTUSubSampleMean) <- c("Set", "numOTU", "sd")
uniOTUSubSampleMean$Day <- gsub(" .*", "", uniOTUSubSampleMean$Set) #regex remove treatment
uniOTUSubSampleMean$Treatment <- gsub("D.* ", "", uniOTUSubSampleMean$Set) #regex day treatment
uniOTUSubSampleMean$Day = factor(uniOTUSubSampleMean$Day, levels = c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D36", "D42"))

View(uniOTUSubSampleMean)

plotMeanOTUBB <- ggplot(uniOTUSubSampleMean, aes(x=Treatment, y=`numOTU`, fill=Treatment), position=position_dodge()) + geom_col( aes(x=Treatment, y=`numOTU`)) + 
  geom_errorbar(aes(ymin=numOTU-(sd), ymax=numOTU+(sd)), width=.5, position=position_dodge(.9)) +
  facet_grid(~Day, scales="free") + theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))  +  ylab("Mean number of distinct OTUs per sample")
plotMeanOTUBB

?geom_boxplot

aov()
ggsave("MeanNumDistinctOTU.tiff", plotMeanOTUBB, width = 15, height = 7, dpi = 500, units =c("in"))
# #D0_IAV_sum$avg <- (as.numeric(D0_IAV_sum$x*100/7))
# uniOTUSubSampleMean$Day = factor(uniOTUSubSampleMean$Day, levels = c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D36", "D42"))
# 
# ggplot(uniOTUSubSampleMean, aes(x=Set, y=numOTU, fill=Day)) + geom_col() + facet_grid(~Day, scales = "free_x") +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + ylab("Total Distinct OTUs on Each Day by Treatment Group")
# ggsave("./TotalDistinctOTUsonEachDaybyTreatmentGroup.tiff", plot=nasalgen.3plot, width = 15, height = 7, dpi = 500, units =c("in"))
# 
# #Mean per group and Day
# uniOTUSubSampleMeanIAV <- uniOTUSubSampleMean %>% filter(Treatment=="IAV")
# uniOTUSubSampleMeanIAV$mean <- (uniOTUSubSampleMeanIAV$numOTU)/7
# View(uniOTUSubSampleMeanIAV)
# 
# uniOTUSubSampleMeanIAV <- uniOTUSubSampleMeanIAV %>% subset(select=-group)
# uniOTUSubSampleMeanIAV <- uniOTUSubSampleMeanIAV %>% distinct()
# 
# plotMeanOTUIAV <- ggplot(uniOTUSubSampleMeanIAV, aes(x=Treatment, y=`mean`, fill=Day)) + geom_col( aes(x=Treatment, y=`mean`, fill=Day)) + facet_grid(~Day, scales = "free_x") +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) 
# 
# #control
# uniOTUSubSampleMeanControl <- uniOTUSubSampleMean %>% filter(Treatment=="Control")
# uniOTUSubSampleMeanControl$mean <- (uniOTUSubSampleMeanControl$numOTU)/10
# 
# 
# uniOTUSubSampleMeanControl <- uniOTUSubSampleMeanControl %>% subset(select=-Pig)
# uniOTUSubSampleMeanControl <- uniOTUSubSampleMeanControl %>% distinct()
# # MeanDistinctOTU <- merge(uniOTUSubSampleMean, uniOTUSubSampleMeanControl, by="Set", ALL.x=TRUE)
# # MeanDistinctOTU <- merge(uniOTUSubSampleMean, uniOTUSubSampleMeanIAV, by="Set", all.x =TRUE)
# plotMeanOTUControl <- ggplot(uniOTUSubSampleMeanControl, aes(x=Treatment, y=`mean`, fill=Day)) + geom_col( aes(x=Treatment, y=`mean`, fill=Day)) + facet_grid(~Day, scales = "free_x") +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) 
# 
# 
# #Combining plots
# otUMeanPlot <- plot_grid(plotMeanOTUIAV + theme(legend.position = "none"), plotMeanOTUControl, labels = "AUTO", rel_widths = c(45.5,54.5), scale=0.9)
# otUMeanPlot