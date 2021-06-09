#########################################
#SRD129 16S - Beta diversity
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
bb.sam <- data.frame(phyloseqbb@sam_data) #Make 'phyloseqbb sam_data' into dataframe. Carry over Phyloseqbb object from 2_phyloseq_BB.R
bb.otu <- data.frame(t(phyloseqbb@otu_table)) #Make 'phyloseqbb otu_table' into dataframe
class(bb.sam) #data.frame
rownames(bb.sam) == row.names(bb.otu) #For rows with sums greater than 1 in 'bb.otu', move rows and their respective sum values into "numOTUs" column in 'bb.sam'
bb.sam$numOTUS <- rowSums(bb.otu > 1)
head(bb.sam)

#NMDS calculation
bb.otu[1:10,1:10]
bb.NMDS <- NMDS_ellipse(bb.sam, bb.otu, grouping_set = 'All')
#Output:
#Result: [1] "Ordination stress: 0.167583915585542"

#Separate meta data and ellipse data to two lists to make NMDS plot
bb.NMDS.2 <- bb.NMDS[[1]]
#'bb.NMDS.2' has meta data + MDS calculations. Select this 1st list of 'bb.NMDS' using double brackets
bb.df_ell.2 <- bb.NMDS[[2]]
#'bb.df_ell.2' is accessing 2nd list from 'bb.NMDS' that has ellipse calculations
#Need two separate lists for making NMDS plot
bb.df_ell.2$group
head(bb.df_ell.2)

#Create "Day" and "Treatment" columns within 'nw.df_ell' for faceting purposes
bb.df_ell.2$Day <- sub(' [A-Za-z]+', '\\1', bb.df_ell.2$group) #Created "Day" column, '\\1' returns the first part of the regular expression D[A-z0-9]+ from 'bb.df_ell.2$group'
bb.df_ell.2$Treatment <- sub('D[A-z0-9]+ ', '\\2', bb.df_ell.2$group) #Create "Treatment" column, '\\2' returns the second part of the sub expression ([A-Za-z]+) from 'bb.df_ell.2$group'
head(bb.df_ell.2)

#Restructure level order for 'nw.metanmds' and 'nw.df_ell'
bb.NMDS.2$Day = factor(bb.NMDS.2$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
bb.df_ell.2$Day = factor(bb.df_ell.2$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
levels(bb.df_ell.2$Day) # [1] "D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"
levels(bb.NMDS.2$Day) # [1] "D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"

#Creating plot from NMDS calculations
bbNMDSplot <- ggplot(data=bb.NMDS.2, aes(x=MDS1, y=MDS2, color=Treatment)) + geom_point() +
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = 0.5) +
  geom_path(data=bb.df_ell.2, aes(x=NMDS1, y=NMDS2, color=Treatment, group=group)) +
  facet_wrap(~Day, scales = 'free') +
  #scale_color_brewer(palette="Dark2") +
  theme_gray(base_size = 10) +
  theme(strip.text.x = element_text(size=15)) +
  labs(caption = 'Ordination stress = 0.17', color="Treatment group")
bbNMDSplot

#Save 'bbNMDSplot' as a .tiff for publication, 500dpi
ggsave("Figure_1.tiff", plot=bbNMDSplot, width = 10, height = 6, dpi = 500, units =c("in"))

#Using pairwise.adonis function
bb.adon <- pairwise.adonis(bb.otu, bb.sam$All, sim.method = 'bray', p.adjust.m = 'bonferroni')
#Run pairwise.adonis on 'bb.otu' OTU table and "All" column of 'bb.sam' dataframe
bb.adon$pairs #List all comparisons in the "pairs" column of 'bb.adon'
goodcomps.bb.adon <- c(grep('D0 [A-Za-z]+ vs D0 [A-Za-z]+', bb.adon$pairs),
                        grep('D1 [A-Za-z]+ vs D1 [A-Za-z]+', bb.adon$pairs),
                        grep('D3 [A-Za-z]+ vs D3 [A-Za-z]+', bb.adon$pairs),
                        grep('D7 [A-Za-z]+ vs D7 [A-Za-z]+', bb.adon$pairs),
                        grep('D10 [A-Za-z]+ vs D10 [A-Za-z]+', bb.adon$pairs),
                        grep('D14 [A-Za-z]+ vs D14 [A-Za-z]+', bb.adon$pairs),
                        grep('D21 [A-Za-z]+ vs D21 [A-Za-z]+', bb.adon$pairs),
                        grep('D28 [A-Za-z]+ vs D28 [A-Za-z]+', bb.adon$pairs),
                        grep('D36 [A-Za-z]+ vs D36 [A-Za-z]+', bb.adon$pairs),
                        grep('D42 [A-Za-z]+ vs D42 [A-Za-z]+', bb.adon$pairs))
# "[A-Za-z]" matches all capital and lowercase letters
# "+" matches a whole word and not just one letter (if you didn't have "+", then it would match by one letter)
# "c" creates the vector, lumps all pairs of specific groups of interest together
# You want to make a vector to combine all the pairwise comparisons you're interested in (same day, different treatment group)
goodcomps.bb.adon
goodcomps.bb.adon.2 <- bb.adon[goodcomps.bb.adon,] #Rename 'goodcomps.bb.adon' vector to 'goodcomps.bb.adon.2'
goodcomps.bb.adon.2
goodcomps.bb.adon.2$p.adjusted <- p.adjust(goodcomps.bb.adon.2$p.value, method = 'fdr') #"p.adjust" function returns a set of p-values adjusted with "fdr" method
goodcomps.bb.adon.2$p.adjusted2 <- round(goodcomps.bb.adon.2$p.adjusted, 3) #Round p-values to 3 decimal points and list in new "p.adjusted2" column
goodcomps.bb.adon.2$p.adjusted2[goodcomps.bb.adon.2$p.adjusted2 > 0.05] <- NA #For all p-values greater than 0.05, replace with "NA"
write.csv(goodcomps.bb.adon.2, file='bb.pairwisecomparisons.csv', row.names=TRUE)

#Alpha diversity: average shannon, invsimpson, numOTUs per "All" subtype
bb.sam$shannon <- diversity(bb.otu) #diversity is vegan function, default index is shannon
#added a shannon index column in nw.meta.nobad
bb.sam$invsimpson <- diversity(bb.otu,index = 'invsimpson')
#use invsimpson since easier to understand than simpson
#don't need to "inverse" simpson values
#simpson: low the number = higher diversity

levels(sample_data(bb.sam)$Day) #"D0"  "D1"  "D10" "D14" "D21" "D28" "D3"  "D36" "D42" "D7"
bb.sam$Day = factor(bb.sam$Day, levels=c("D0", "D1", "D3", "D7", "D10", "D14", "D21", "D28", "D36", "D42"))
levels(sample_data(bb.sam)$Day) #"D0"  "D1"  "D3"  "D7"  "D10" "D14" "D21" "D28" "D36" "D42"   
head(bb.sam)

#average shannon, invsimpson, numOTUs per "All" subtype
bb.sam.average.shannon.invsimpson.numOTUs <- aggregate(bb.sam[, 6:8], list(bb.sam$All), mean)
print(bb.sam.average.shannon.invsimpson.numOTUs)

#         Group.1   numOTUS  shannon invsimpson
#1        D0 BB  43.00000 2.023254   4.419315
# 2   D0 Control  70.20000 2.553821   6.328268
# 3        D1 BB  47.55556 2.268927   5.997870
# 4   D1 Control  51.00000 2.253296   5.126313
# 5       D10 BB  71.60000 2.646794   7.370097
# 6  D10 Control  62.50000 2.442391   5.105919
# 7       D14 BB  73.30000 2.496178   4.141793
# 8  D14 Control  81.80000 2.855901   7.141871
# 9       D21 BB  89.10000 3.096551   7.973830
# 10 D21 Control  79.60000 2.933431   7.151400
# 11      D28 BB  96.50000 3.312364  11.275079
# 12 D28 Control 104.25000 3.464386  13.549240
# 13       D3 BB  86.88889 2.916247   7.478663
# 14  D3 Control  82.90000 2.889451   7.220153
# 15      D36 BB  97.00000 3.416072  13.706073
# 16 D36 Control  91.30000 3.147094  12.772148
# 17      D42 BB 129.90000 3.906807  28.307299
# 18 D42 Control 106.30000 3.488465  17.076775
# 19       D7 BB  63.62500 2.362948   5.267264
# 20  D7 Control  59.70000 2.393280   5.465053

bb.pairwise.wilcox.shannon.test <- pairwise.wilcox.test(bb.sam$shannon, bb.sam$All, p.adjust.method = 'none') #Calculate pairwise comparisons by "All" column of the shannon indices in "Shannon" column
print(bb.pairwise.wilcox.shannon.test) #Look at the results of 'nw2.pairwise.wilcox.shannon.test'
bb.pairwise.wilcox.invsimpson.test <- pairwise.wilcox.test(bb.sam$invsimpson, bb.sam$All, p.adjust.method = 'none')
print(bb.pairwise.wilcox.invsimpson.test)

#box and whisker plot: shannon
bb.shan.bw <- ggplot(data = bb.sam, aes(x=Treatment, y=shannon, group=All, fill=Treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~Day, scales = 'free') +
  scale_y_continuous(name = "Shannon diversity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text.x = element_text(size=10),
        axis.title.y = element_text(size=15))
#free allows each plot to customize scale to the specific data set (no forced scale applied to all plots)
#position = position_dodge2(preserve = 'total') --> fixes the ggplot box width;
#prevents narrow boxes from forming in plot (they become wider)
bb.shan.bw

#box and whisker plot: inverse simpson
bb.invsimp.bw <- ggplot(data = bb.sam, aes(x=Treatment, y=invsimpson, group=All, fill=Treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~Day, scales = 'free') +
  scale_y_continuous(name = "Inverse Simpson diversity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text.x = element_text(size=10),
        axis.title.y = element_text(size=15))
bb.invsimp.bw

#combine shannon and inverse simpson box and whisker plots
bb.combalpha <- plot_grid(bb.shan.bw + theme(legend.position = "none"), bb.invsimp.bw, labels = "AUTO")
bb.combalpha

#Save 'bb.combalpha' as a .tiff for publication, 500dpi
ggsave("SRD129_BBControl_ShannonInverseSimpsonCombined.tiff", plot=bb.combalpha, width = 10, height = 5, dpi = 500, units =c("in"))
