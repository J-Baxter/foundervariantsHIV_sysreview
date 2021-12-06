###################################################################################################
###################################################################################################
#Figures describing covariate distirbution accompanying Inferring the multiplicity of founder variants
# initiating HIV-1 infection: a systematic review and IPD meta-analysis
# 1. Panel of bar plots (and one histogram) plotting frequency against covariate, stacking bars by 
# binary classification of multiplicity
###################################################################################################
###################################################################################################

# RUN FROM HERE #
# Dependencies
source('./load_packages.R')
source('./scripts/generalpurpose_funcs.R')


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
if (!dir.exists('data')){
  Retrieve('data.zip')
}else{
  Sys.sleep(0.2)
}

df <- read.csv("./data/meta_analysis_data.csv",
               na.strings = "NA",
               stringsAsFactors = T) %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene','sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  droplevels()

# Set reference levels for meta regression
# HSX:MTF, haplotype (highlighter), unknown seropositivity, B, whole genome
baseline.covar <- c("reported.exposure_", "grouped.method_", "grouped.subtype_","sequencing.gene_", "sampling.delay_",'alignment.bin_', 'sequencing.method_')
baseline.level <- c("HSX:MTF", "haplotype", "B" , "whole.genome" , "<21", 'NFLG', 'sanger_SGA')

df <- SetBaseline(df, baseline.covar, baseline.level)

###################################################################################################
###################################################################################################
# Basic plots for figure 2 A-H (NB figure 1 is the PRISMA flowchart)

# Set colour palettes 
mycols_founder <- RColorBrewer::brewer.pal(name = 'RdBu', n = 8)[c(2,7)] #c("#E64B35FF", "#4DBBD5FF")
nb.cols <- 12
mycols_method <- colorRampPalette(brewer.pal(10, "RdBu"))(nb.cols)

# Set Labels
labs <- c('Multiple','Single')

###################################################################################################
# Fig 2A: Detailed barplot displaying route of transmission (with direction subcategories)
colnames <- c('reported.exposure' , 'sub.exposure' , 'frequency')
exposures_df <- stacked_categories(df$reported.exposure_, colnames)
exposures_df$sub.exposure <- factor(exposures_df$sub.exposure, levels = c('MTF', 'FTM', 'nodirection',
                                                                          'MSM',
                                                                          'PreP', 'IntraP', 'PostP', 'notiming',
                                                                          'PWID'))

fig2_a <- ggplot(exposures_df, aes(x = reported.exposure, y = frequency))+
  geom_bar(stat = 'identity' , aes(fill = sub.exposure, position = 'stack'))+
  scale_fill_manual(values= mycols_method ,labels = c('HSX: MTF', 'HSX: FTM', 'HSX: undisclosed',
                                                      'MSM',
                                                      'MTC: Pre-Partum', 'MTC: Intrapartum', 'MTC: Post-Partum', 'MTC: undisclosed',
                                                      'PWID'))+
  scale_y_continuous(limits = c(0,1000), expand = c(0,0)) +
  theme_classic() +
  xlab('Risk Group') +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Number of Participants') +
  labs(fill = "Reported Exposure", colour = "Reported Exposure")


###################################################################################################
# Fig 2B: Route of transmission (riskgroup)
fig2_b <- ggplot(df, aes(x = riskgroup_))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_)), y = (..count..)/sum(..count..), colour = forcats::fct_rev(factor(multiple.founders_))))+
  scale_color_manual(values = mycols_founder, labels = labs)+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_discrete(labels = LabelX(df$riskgroup_)) +
  theme_classic()+
  xlab('Risk Group')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Proportion of Participants')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity") 


###################################################################################################
# Fig 2C: Method of quantification
fig2_c <- ggplot(df, aes(x = grouped.method_))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_)), y = (..count..)/sum(..count..)))+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
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
  ylab('Proportion of Participants')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")


###################################################################################################
# Fig 2D: Subtype
fig2_d <- ggplot(df, aes(x = grouped.subtype_))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_)), y = (..count..)/sum(..count..)))+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme_classic()+
  xlab('Subtype')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Proportion of Participants')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")


###################################################################################################
# Fig 2E: Sampling Delay (inferred from feibig/seropositivty/timing of sampling relative to infection)
fig2_e <- ggplot(df, aes(x = sampling.delay_))+ 
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_)), y = (..count..)/sum(..count..)))+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_discrete(labels = c('<21' = '<21 Days', '>21' = '>21 Days', 'Unknown' = 'Unknown') %>% str_to_title()) +
  theme_classic()+
  xlab('Sampling Delay')+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Proportion of Participants')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")


###################################################################################################
# Fig 2F: Number of consensus sequences analysed
numseqs_df <- GetNumSeqs(df)

fig2_f <- ggplot(numseqs_df , aes(sequencing.number_))+
  geom_histogram(aes(fill = forcats::fct_rev(factor(multiple.founders_)), y = (..count..)/sum(..count..)), binwidth = 3)+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,90), expand = c(0,0))+
  theme_classic()+
  xlab('Number of Consensus Sequences')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Proportion of Participants')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")



###################################################################################################
# Fig 2G: Genomic region analysed

fig2_g <- ggplot(df, aes(x = sequencing.gene_))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_)), y = (..count..)/sum(..count..)))+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_discrete(labels = c('env' = 'Env',
                              'gag'= 'Gag',
                              'pol'='Pol',
                              'whole.genome' = 'NFLG')) +
  theme_classic()+
  xlab('Genomic Region Analysed')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Proportion of Participants')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")


###################################################################################################
# Fig 2H =  Alignment Length

fig2_h <- ggplot(df, aes(alignment.length_))+
  geom_histogram(aes(fill = forcats::fct_rev(factor(multiple.founders_)), y = (..count..)/sum(..count..)), binwidth = 100)+
  scale_fill_manual(values = mycols_founder, labels = labs)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,9500), expand = c(0,0)) +
  theme_classic()+
  xlab('Alignment Length')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Proportion of Participants') +
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity") 


###################################################################################################
# Combine basic covariate plots into figure (with labels and axis)
fig2_top <- cowplot::plot_grid(fig2_a + theme(legend.position= c(0.85,0.835)),
                               fig2_b + theme(legend.position= c(0.85,0.95)),
                               nrow = 1,
                               align = "hv",
                               axis = "bt" ,
                               labels = "AUTO")

fig2_bottom <- cowplot::plot_grid(fig2_c + theme(legend.position="none"),
                                  fig2_d + theme(legend.position="none", axis.title.x = element_text(margin = margin(t = 13))), 
                                  fig2_e + theme(legend.position="none", axis.title.x = element_text(margin = margin(t = 20))),
                                  fig2_f + theme(legend.position="none", axis.title.x = element_text(margin = margin(t = 12))),
                                  fig2_g + theme(legend.position="none"),
                                  fig2_h + theme(legend.position="none", axis.title.x = element_text(margin = margin(t = 6))),
                                  ncol = 3, 
                                  nrow = 2,
                                  align = "hv", 
                                  axis = "bt",
                                  labels = c('C', 'D', 'E', 'F', 'G', 'H')) 

fig_2 <- cowplot::plot_grid(fig2_top, fig2_bottom, nrow = 2)

# Print to file
setEPS()
postscript("./results/figure2.eps", width = 10, height = 16)
fig_2
dev.off()

###################################################################################################
# Presentation plots
legend <- get_legend(p1.method )
prow <- cowplot::plot_grid(p1.method + theme(legend.position="none"),
                           p1.subtype + theme(legend.position="none"), 
                           p1.seropos + theme(legend.position="none"),
                           p1.numseq + theme(legend.position="none"),
                           p1.sg + theme(legend.position="none"),
                           p1.alignment + theme(legend.position="none"),
                           ncol = 3 , align = "hv", axis = "bt",
                           labels = "AUTO", rel_widths =c(1, 1))

cowplot::plot_grid(prow, legend, ncol = 2, rel_widths  = c(1,0.1))




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
# Panel S4: Geographic and Timewise structere of transmission

# Deprecated #
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
  
  theme_minimal(base_size = 10) +
  
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
###################################################################################################

# Year of Publication ~ Frequency of Individuals, stacked by risk group
figureS4_a <- ggplot(df, aes(year_))+
  geom_histogram(aes(fill = riskgroup_, 
                     colour = riskgroup_), 
                 bins = 8) +
  scale_color_manual(values = mycols_method[c(2,4,6,8,10,12)]) +
  scale_fill_manual(values = mycols_method[c(2,4,6,8,10,12)]) +
  scale_y_continuous(limits = c(0,500), 
                     expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic(base_size = 10)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  
  labs(fill = 'Risk Group', 
       x = 'Year of Publication', 
       y = 'Proportion of Participants')+
  guides(color = FALSE)+
  theme(legend.position = 'bottom')

# Year of Publication ~ Frequency of Individuals, stacked by method of enumeration
figureS4_b <- ggplot(df, aes(year_))+
  geom_histogram(aes(fill = grouped.method_, 
                     colour = grouped.method_), 
                 bins = 8) +
  scale_color_manual(values = mycols_method[c(2,4,6,8,10,12)]) +
  scale_fill_manual(values = mycols_method[c(2,4,6,8,10,12)],
                    labels = c('Haplotype', 
                               'Distance', 
                               'Model',
                               'Molecular',
                               'Phylogenetic: Recipient Only',
                               'Phylogenetic: Source & Recipient')) +
  scale_y_continuous(limits = c(0,500), 
                     expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic(base_size = 10)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  
  labs(fill = 'Method', 
       x = 'Year of Publication', 
       y = 'Proportion of Participants')+
  guides(color = FALSE)+
  theme(legend.position = 'bottom') 

# Year of Publication ~ Frequency of Individuals, stacked by sequencing technology
df$sequencing.method_<- factor(df$sequencing.method_, levels = c('sanger_SGA', 
                                                                 '2G:illumina_miseq', 
                                                                 '2G:iontorrent',
                                                                 '2G:roche_454',
                                                                 '3G:PacBio_hifi',
                                                                 'sanger', 
                                                                 'sanger_precSGA',
                                                                 'unknown'))


figureS4_c <- ggplot(df, aes(year_))+
  geom_histogram(aes(fill = sequencing.method_, 
                     colour = sequencing.method_), 
                 bins = 8) +
  scale_color_manual(values = mycols_method[c(2,4,6,8,10,12,1,11)]) +
  scale_fill_manual(values = mycols_method[c(2,4,6,8,10,12,1,11)],
                    labels = c('Sanger with SGA', 
                               'Illumina MiSeq', 
                               'ONT IonTorrent',
                               'Roche 454',
                               'PacBio HiFi',
                               'Sanger Only', 
                               'Sanger with SGA Precursor',
                               'Unknown')) +
  scale_y_continuous(limits = c(0,500), 
                     expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic(base_size = 10)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(fill = 'Sequencing Technology', 
       x = 'Year of Publication', 
       y = 'Proportion of Participants')+
  guides(color = FALSE, fill = guide_legend(nrow = 3, byrow= TRUE))+
  theme(legend.position = 'bottom') 

figureS4 <- cowplot::plot_grid(figureS4_a, 
                               figureS4_b + theme(axis.title.y = element_blank()), 
                               figureS4_c + theme(axis.title.y = element_blank()), 
                               ncol = 3, labels = 'AUTO', align = 'hv')

# Print to file
setEPS()
postscript("./results/figureS4.eps", width = 18, height = 8)
figureS4
dev.off()

###################################################################################################
###################################################################################################
# Delay vs Method ~ deprecated
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