###################################################################################################
###################################################################################################
#Figures describing covariate distirbution accompanying Inferring the multiplicity of founder variants
# initiating HIV-1 infection: a systematic review and IPD meta-analysis
# 1. Panel of bar plots (and one histogram) plotting frequency against covariate, stacking bars by 
# binary classification of multiplicity
###################################################################################################
###################################################################################################

# Dependencies
library(ggplot2)
library(ggsci)
library(dplyr)
library(forcats)
library(tidyr)
library(RColorBrewer)
library(cowplot)
source('generalpurpose_funcs.R')


#function for stacking categories and calculating summary frequencies
stacked_categories <- function(x, catnames){
  characters <- as.character(x)
  split_list <- strsplit(characters , '[:]')
  split_df <- do.call(rbind.data.frame, split_list)
  
  freq <- split_df %>%
    group_by(split_df[,1] , split_df[,2]) %>%
    summarise(frequency = n())
  
  colnames(freq) <- catnames
  return(freq)
}


# Extracting and formatting sequence numbers from main dataframe
GetNumSeqs <- function(df){
  # drop NA and NGS first 
  numseqs <- df[!df$sequencing.number %in% c('NGS', NA),]
  numseqs$sequencing.number <- as.numeric(numseqs_df$sequencing.number)
  
  return(numseqs)
}


###################################################################################################
###################################################################################################
# Import data
setwd("./data")
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% formatDF()


###################################################################################################
###################################################################################################
# Basic plots for panel 1 (four bar plots to panel displaying frequencies of grouped method, 
# seropositivity, subtype and number of sequences)

# Set colour palettes 
mycols_founder <- RColorBrewer::brewer.pal(name = 'RdBu', n = 8)[c(2,7)]
nb.cols <- 12
mycols_method <- colorRampPalette(brewer.pal(10, "RdBu"))(nb.cols)


# Set Labels
labs <- c('Multiple','Sinlge')


###################################################################################################
# 1.1 Method
p1.method <- ggplot(df, aes(grouped.method))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders)), colour = forcats::fct_rev(factor(multiple.founders))))+
  scale_color_manual(values = mycols_founder, labels = labs)+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  ylim(0,1000)+
  theme_classic()+
  xlab('Grouped Method')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")


###################################################################################################
# 1.2 Seroconversion (poss stack infant and NA together)
p1.seropos <- ggplot(df, aes(participant.seropositivity))+ 
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders)), colour = forcats::fct_rev(factor(multiple.founders))))+
  scale_color_manual(values = mycols_founder, labels = labs)+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  ylim(0,1000)+
  theme_classic()+
  xlab('Seropositivity')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")


###################################################################################################
# 1.3 Subtype
p1.subtype <- ggplot(df, aes(grouped.subtype))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders)), colour = forcats::fct_rev(factor(multiple.founders))))+
  scale_color_manual(values = mycols_founder, labels = labs)+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  ylim(0,1000)+
  theme_classic()+
  xlab('Subtype')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")


###################################################################################################
# 1.4 Number of consensus sequences analysed
numseqs_df <- GetNumSeqs(df)

p1.numseq <- ggplot(numseqs_df , aes(sequencing.number))+
  geom_histogram(binwidth=5, (aes(fill = forcats::fct_rev(factor(multiple.founders)), colour = forcats::fct_rev(factor(multiple.founders)))))+
  scale_color_manual(values = mycols_founder, labels = labs)+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  ylim(0,1000)+
  theme_classic()+
  xlab('Number of Consensus Genomes Analysed')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency') +
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")


###################################################################################################
# 1.5 Route of transmission (reported.exposure)

p1.exposure <- ggplot(df, aes(riskgroup))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders)), colour = forcats::fct_rev(factor(multiple.founders))))+
  scale_color_manual(values = mycols_founder, labels = labs)+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  ylim(0,1000)+
  theme_classic()+
  xlab('Reported Exposure')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity") 


###################################################################################################
# 1.6 Detailed barplot displaying route of transmission (with direction subcategories)
names <- c('reported.exposure' , 'sub.exposure' , 'frequency')
exposures_df <- stacked_categories(df$reported.exposure, names)

p1.6 <- ggplot(exposures_df, aes(x = reported.exposure , y = frequency))+
  geom_bar(stat = 'identity' , aes(fill = sub.exposure, colour = sub.exposure), position = 'stack')+
  scale_color_manual(values = mycols_method)+
  scale_fill_manual(values = mycols_method)+
  theme_classic()+
  xlab('Exposure')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency') + labs(fill = "Reported Exposure", colour = "Reported Exposure")


###################################################################################################
# Combine basic covariate plots into figure (with labels and axis)
plot_grid(p1.method + theme(legend.position="none"),
          p1.subtype + theme(legend.position="none"), 
          p1.seropos + theme(legend.position="none"),
          p1.numseq + theme(legend.position="none"),
          p1.exposure + theme(legend.position="none"),
          p1.6,
          ncol = 3 , nrow = 3,align = "hv", axis = "bt" , labels = "AUTO") 

# Print to file
tiff("covar_barplot.tiff" , width = 14 , height = 20)
figure1
dev.off()

###################################################################################################
###################################################################################################
# END #
###################################################################################################
###################################################################################################