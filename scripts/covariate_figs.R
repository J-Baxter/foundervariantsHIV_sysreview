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
library(stringr)
source('~/foundervariantsHIV_sysreview/scripts/generalpurpose_funcs.R')


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
GetNumSeqs <- function(data){
  # drop NA and NGS first 
  numseqs <- data[!data$sequencing.number_ %in% c('NGS', NA),]
  numseqs$sequencing.number_ <- as.numeric(numseqs$sequencing.number_)
  
  return(numseqs)
}


# Format categorical X axis tick lables
LabelX <- function(xvar, caps){
  if (class(xvar) == 'factor'){
    labs <- gsub("[.]" , " ", levels(xvar)) %>%
      str_wrap(.,width = 10)
  } else{
    warning('variable is not discrete')
  }
  
  return(labs)
}
###################################################################################################
###################################################################################################
# Import data (as for meta regression)
df <- read.csv("./data/data_master_11121.csv", na.strings = "NA") %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene','sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  droplevels()


###################################################################################################
###################################################################################################
# Basic plots for panel 1 (four bar plots to panel displaying frequencies of grouped method, 
# seropositivity, subtype and number of sequences)

# Set colour palettes 
mycols_founder <- RColorBrewer::brewer.pal(name = 'RdBu', n = 8)[c(2,7)] #c("#E64B35FF", "#4DBBD5FF")

nb.cols <- 12
mycols_method <- colorRampPalette(brewer.pal(10, "RdBu"))(nb.cols)


# Set Labels
labs <- c('Multiple','Single')


###################################################################################################
# 1.1 Method
p1.method <- ggplot(df, aes(grouped.method_))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_))))+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1350), expand = c(0,0)) +
  scale_x_discrete(labels = c('distance' = 'Distance',
                              'haplotype' = 'Haplotype',
                              'model' = 'Model', 
                              'molecular' = 'Molecular',
                              'phylogenetic:r' = 'Phylogenetic: recipient only',
                              'phylogenetic:s&R '= 'Phylogenetic: paired') %>%
                     str_wrap(width = 17)) +
  theme_classic()+
  xlab('Method of Quantification')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")


###################################################################################################
# 1.2 Sampling Delay (inferred from feibig/seropositivty/timing of sampling relative to infection)
p1.seropos <- ggplot(df, aes(sampling.delay_))+ 
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_))))+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1350), expand = c(0,0)) +
  scale_x_discrete(labels = c('<21' = '<21 Days', '>21' = '>21 Days', 'Unknown' = 'Unknown') %>% str_to_title()) +
  theme_classic()+
  xlab('Sampling Delay')+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")


###################################################################################################
# 1.3 Subtype
p1.subtype <- ggplot(df, aes(grouped.subtype_))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_))))+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1350), expand = c(0,0)) +
  theme_classic()+
  xlab('Subtype')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")


###################################################################################################
# 1.4 Number of consensus sequences analysed
numseqs_df <- GetNumSeqs(df)

p1.numseq <- ggplot(numseqs_df , aes(sequencing.number_))+
  geom_histogram(binwidth=5, (aes(fill = forcats::fct_rev(factor(multiple.founders_)))))+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1350), expand = c(0,0)) +
  theme_classic()+
  xlab('Number of Consensus Sequences')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency') +
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")


###################################################################################################
# 1.5 Route of transmission (riskgroup)

p1.exposure <- ggplot(df, aes(riskgroup_))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_)), colour = forcats::fct_rev(factor(multiple.founders_))))+
  scale_color_manual(values = mycols_founder, labels = labs)+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1000), expand = c(0,0)) +
  scale_x_discrete(labels = LabelX(df$riskgroup_)) +
  theme_classic()+
  xlab('Risk Group')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity") 


###################################################################################################
# 1.6 Detailed barplot displaying route of transmission (with direction subcategories)
colnames <- c('reported.exposure' , 'sub.exposure' , 'frequency')
exposures_df <- stacked_categories(df$reported.exposure_, colnames)
exposures_df$sub.exposure <- factor(exposures_df$sub.exposure, levels = c('MTF', 'FTM', 'nodirection',
                                                                          'MSM',
                                                                          'PreP', 'IntraP', 'PostP', 'notiming',
                                                                          'PWID'))

p1.6 <- ggplot(exposures_df, aes(x = reported.exposure , y = frequency))+
  geom_bar(stat = 'identity' , aes(fill = sub.exposure), position = 'stack')+
  scale_fill_manual(values= mycols_method ,labels = c('HSX: MTF', 'HSX: FTM', 'HSX: undisclosed',
                                                      'MSM',
                                                      'MTC: Pre-Partum', 'MTC: Intrapartum', 'MTC: Post-Partum', 'MTC: undisclosed',
                                                      'PWID'))+
  scale_y_continuous(limits = c(0,1000), expand = c(0,0)) +
  theme_classic()+
  xlab('Risk Group')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency') + labs(fill = "Reported Exposure", colour = "Reported Exposure")

###################################################################################################
# 1.7 Sequencing Gene

p1.sg <- ggplot(df, aes(sequencing.gene_))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_))))+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1350), expand = c(0,0)) +
  scale_x_discrete(labels = c('env' = 'Env',
                              'gag'= 'Gag',
                              'pol'='Pol',
                              'whole.genome' = 'NFLG')) +
  theme_classic()+
  xlab('Gene Analysed')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")



###################################################################################################
# 1.8 Alignment Length

p1.alignment <- ggplot(df, aes(alignment.length_))+
  geom_histogram(aes(fill = forcats::fct_rev(factor(multiple.founders_))))+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1350), expand = c(0,0)) +
  theme_classic()+
  xlab('Alignment Length')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity") 

###################################################################################################
# Combine basic covariate plots into figure (with labels and axis)
grid1 <- cowplot::plot_grid( p1.6+theme(legend.position= c(0.85,0.7)),
                             p1.exposure + theme(legend.position= c(0.85,0.87)),
                             nrow = 2 ,align = "hv", axis = "bt" , labels = "AUTO", rel_heights = c(2,1))

grid2 <- cowplot::plot_grid(p1.method + theme(legend.position="none"),
                            p1.subtype + theme(legend.position="none"), 
                            p1.seropos + theme(legend.position="none"),
                            p1.numseq + theme(legend.position="none"),
                            p1.sg + theme(legend.position="none"),
                            p1.alignment + theme(legend.position="none"),
                            ncol = 3 , nrow = 2,align = "hv", axis = "bt" , labels = c("C", "D" , "E" , "F" , "G" , "H")) 

cowplot::plot_grid(grid1, grid2, nrow = 2)

legend <- get_legend(p1.method )
# Print to file
jpeg("./results/panel1a.jpeg" ,width = 6000, height = 4000, res = 380 ,units = "px", pointsize = 12)

prow <- cowplot::plot_grid(p1.method + theme(legend.position="none"),
                           p1.subtype + theme(legend.position="none"), 
                           p1.seropos + theme(legend.position="none"),
                           p1.numseq + theme(legend.position="none"),
                           p1.sg + theme(legend.position="none"),
                           p1.alignment + theme(legend.position="none"),
                           ncol = 3 , align = "hv", axis = "bt",
                           labels = "AUTO", rel_widths =c(1, 1))

cowplot::plot_grid(prow, legend, ncol = 2, rel_widths  = c(1,0.1))
dev.off()

library(gridExtra)
grid.arrange(p1.exposure + theme(legend.position= c(0.85,0.87)), 
             p1.6+theme(legend.position= c(0.85,0.7)),
             p1.method + theme(legend.position="none"),
             p1.subtype + theme(legend.position="none"), 
             p1.seropos + theme(legend.position="none"),
             p1.numseq + theme(legend.position="none"),
             p1.sg + theme(legend.position="none"),
             p1.alignment + theme(legend.position="none"), ncol = 3, layout_matrix = cbind(c(2,2,1), c(3,4,5), c(6,7,8)))
###################################################################################################
###################################################################################################
# Geographic and Timewise structere of transmission

region_exp <- df %>%
  group_by(reported.exposure_, grouped.subtype_) %>%
  summarise(subjects = n(), multiplefounders = sum(multiple.founders_)) %>%
  filter(grouped.subtype_ != '') %>%
  droplevels()

p2.1 <- ggplot(region_exp,aes(x = reported.exposure_,
                              y = grouped.subtype_,
                              size = subjects, 
                              color = (multiplefounders/subjects))) +
  geom_point() + 
  
  theme_minimal(base_size = 10,
                base_family = 'sans') +
  
  scale_color_distiller(palette = 'RdBu') +
  
  theme(axis.line = element_blank(),
        legend.position = 'bottom',
        legend.box = 'vertical') +
  
  scale_size(range = c(2,15)) + 
  
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  
  labs(size = 'Number of Patients', 
       color = 'Frequency of Founder Variant Multiplicity', 
       x = 'Reported Exposure', 
       y = 'Virus Subtype')

p2.2 <- ggplot(df, aes(year_))+
  geom_histogram(aes(fill = riskgroup_, 
                     colour = riskgroup_), 
                 bins = 8) +
  scale_color_manual(values = mycols_method[c(2,4,6,8,10,12)]) +
  scale_fill_manual(values = mycols_method[c(2,4,6,8,10,12)]) +
  scale_y_continuous(limits = c(0,500), 
                     expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic(base_size = 10,
                base_family = 'sans')+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  
  labs(fill = 'Risk Group', 
       x = 'Year of Publication', 
       y = 'Frequency')+
  guides(color = FALSE)+
  theme(legend.position = 'bottom') 

p2.3 <- ggplot(df, aes(year_))+
  geom_histogram(aes(fill = grouped.method_, 
                     colour = grouped.method_), 
                 bins = 8) +
  scale_color_manual(values = mycols_method[c(2,4,6,8,10,12)]) +
  scale_fill_manual(values = mycols_method[c(2,4,6,8,10,12)]) +
  scale_y_continuous(limits = c(0,500), 
                     expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic(base_size = 10,
                base_family = 'sans')+
  theme(axis.text.x=element_text(angle=45, hjust=1), size = 10)+
  
  labs(fill = 'Method', 
       x = 'Year of Publication', 
       y = 'Frequency')+
  guides(color = FALSE)+
  theme(legend.position = 'bottom') 

grid3 <- cowplot::plot_grid( p2.2,p2.3, ncol = 2, labels = 'AUTO', align = 'hv')

jpeg("./results/subgroup_transmission.jpeg" ,width = 6000, height = 4000, res = 380 ,units = "px", pointsize = 12)
grid3 
dev.off()


###################################################################################################
###################################################################################################
# Delay vs Method
delay_method <- df %>%
  group_by(sampling.delay_, grouped.method_) %>%
  summarise(subjects = n(), multiplefounders = sum(multiple.founders_)) %>%
  droplevels()


p3.1 <- ggplot(delay_method, aes(x = sampling.delay_,
                                 y = grouped.method_,
                                 size = subjects, 
                                 color = (multiplefounders/subjects))) +
  geom_point() + 
  
  theme_minimal(base_size = 10,
                base_family = 'sans') +
  
  scale_color_viridis_c() +
  
  theme(axis.line = element_blank(),
        legend.position = 'bottom',
        legend.box = 'vertical') +
  
  scale_size(range = c(1,15)) + 
  
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  
  labs(size = 'Number of Patients', 
       color = 'Frequency of Founder Variant Multiplicity', 
       x = 'Sampling Delay', 
       y = 'Method')

delay_exp <- df %>%
  group_by(sampling.delay_, reported.exposure_) %>%
  summarise(subjects = n(), multiplefounders = sum(multiple.founders_)) %>%
  droplevels()


p3.2 <- ggplot(delay_exp, aes(x = sampling.delay_,
                              y = reported.exposure_,
                              size = subjects, 
                              color = (multiplefounders/subjects))) +
  geom_point() + 
  
  theme_minimal(base_size = 10,
                base_family = 'sans') +
  
  scale_color_viridis_c() +
  
  theme(axis.line = element_blank(),
        legend.position = 'bottom',
        legend.box = 'vertical') +
  
  scale_size(range = c(1,15)) + 
  
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  
  labs(size = 'Number of Patients', 
       color = 'Frequency of Founder Variant Multiplicity', 
       x = 'Sampling Delay', 
       y = 'Reported Exposure')
###################################################################################################
###################################################################################################
# END #
###################################################################################################
###################################################################################################