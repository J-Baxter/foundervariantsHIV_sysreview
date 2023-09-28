# Thesis Plots (Chapter two of PhD Thesis)

####################################### Dependencies #######################################
library(tidyverse)
library(ggplot2)
require(extrafont)
source('./scripts/generalpurpose_funcs.R')

#https://www.fontsquirrel.com/fonts/latin-modern-sans
font_import(pattern = "lmsans10*") 
loadfonts()

my_theme <- theme_classic(base_family = "LM Sans 10")+
  theme(
    text = element_text(size=16),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = 'none', 
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank()
  )


#function for stacking categories and calculating summary frequencies
stacked_categories <- function(x, catnames){
  characters <- as.character(x)
  split_list <- strsplit(characters , '[:]')
  split_df <- do.call(rbind.data.frame, split_list)
  
  freq <- split_df %>%
    group_by(split_df[,1] , split_df[,2]) %>%
    dplyr::summarise(frequency = n()) 
  
  
  
  colnames(freq) <- catnames
  
  return(freq)
}


# Extracting and formatting sequence numbers from main dataframe
GetNumSeqs <- function(data){
  # drop NA and NGS first 
  out <- data  %>% 
    filter(!grepl('unknown', sequencing.number_)) %>% 
    mutate(sequencing.number_ = unlist(sequencing.number_) %>% 
             unname() %>% 
             as.character() %>% 
             as.numeric()) %>% 
    filter(!is.na(sequencing.number_))
  
  
  return(out)
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


####################################### Import Data #######################################
df <- read.csv("./data/meta_analysis_data.csv",
               na.strings = "NA",
               stringsAsFactors = TRUE) %>%
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene','sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  droplevels()

# Set reference levels for meta regression
# HSX:MTF, haplotype (highlighter), unknown seropositivity, B, whole genome
baseline.covar <- c("reported.exposure_", "grouped.method_", "grouped.subtype_","sequencing.gene_", "sampling.delay_",'alignment.bin_', 'sequencing.method_')
baseline.level <- c("HSX:MTF", "haplotype", "B" , "env" , "<21", 'NFLG', 'sanger_SGA')

df <- SetBaseline(df, baseline.covar, baseline.level)


####################################### Figure 1 #######################################
# Figure 1 will be the PRISMA flow-chart

####################################### Figure 2 #######################################
labs <- c('Multiple','Single')

colnames <- c('reported.exposure' , 'sub.exposure' , 'frequency')
exposures_df <- stacked_categories(df$reported.exposure_, colnames)
exposures_df$sub.exposure <- factor(exposures_df$sub.exposure, levels = c('MTF', 'FTM', 'nodirection',
                                                                          'MSM',
                                                                          'PreP', 'IntraP', 'PostP', 'notiming',
                                                                          'PWID'))

# Fig 2A
plt_2a <- ggplot(exposures_df,
                 aes(x = reported.exposure, y = frequency))+
  geom_bar(stat = 'identity',
           aes(fill = sub.exposure), 
           position = 'stack')+
  scale_fill_brewer("Reported Exposure",
                    palette = 'YlGn',
                    labels = c('MTF',
                               'FTM', 
                               'HSX: Unknown',
                               'MSM',
                               'MTC: Pre-Partum', 
                               'MTC: Intrapartum', 
                               'MTC: Post-Partum', 
                               'MTC: Unknown',
                               'PWID'))+
  scale_y_continuous('Number of Participants', 
                     limits = c(0,1000), 
                     expand = c(0,0)) +
  scale_x_discrete('Risk Group',
                   expand = c(0,0)) +
  my_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1))


# Fig 2B: Route of transmission (riskgroup)
plt_2b <- ggplot(df, aes(x = reported.exposure_))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_)),
               y = (..count..)/sum(..count..)))+
  scale_fill_brewer("Founder Multiplicity",
                    palette = 'YlGn',
                    labels = labs) +
  scale_y_continuous('Proportion of Participants',
                     limits = c(0,1),
                     expand = c(0,0)) +
  scale_x_discrete('Risk Group',
                   labels = c('MTF', 
                              'FTM',
                              'HSX: Undisclosed',
                              'MSM',
                              'MTC: Intrapartum', 
                              'MTC: Undisclosed', 
                              'MTC: Post-Partum',
                              'MTC: Pre-Partum', 
                              'PWID'),
    expand = c(0,0)) +
  my_theme+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.y = element_blank())

  
# Fig 2C: Method of quantification
plt_2c <- ggplot(df, aes(x = grouped.method_))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_)),
               y = (..count..)/sum(..count..)))+
  scale_fill_brewer("Founder Multiplicity",
                    palette = 'YlGn',
                    labels = labs) +
  scale_y_continuous('Proportion of Participants',
                     limits = c(0,1),
                     expand = c(0,0)) +
  scale_x_discrete('Method of Quantification',
                   labels = c('distance' = 'Distance',
                              'haplotype' = 'Haplotype',
                              'model' = 'Model', 
                              'molecular' = 'Molecular',
                              'phylogenetic:r' = 'Phylogenetic: Recipient Only',
                              'phylogenetic:s&R '= 'Phylogenetic: Paired') %>%
                     str_wrap(width = 17),
                   expand = c(0,0)) +
  my_theme+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.y = element_blank())
  

# Fig 2D: Subtype
subtype_order <- c('A','B', 'C', 'D' , 'CRF01_AE' , 'CRF02_AG' , 'recombinant', 'unknown')
df$grouped.subtype_ <- factor(df$grouped.subtype_, levels = subtype_order)
table(df$grouped.subtype_)

plt_2d <- ggplot(df, aes(x = grouped.subtype_))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_)), 
               y = (..count..)/sum(..count..)))+
  scale_fill_brewer("Founder Multiplicity",
                    palette = 'YlGn',
                    labels = labs) +
  scale_y_continuous('Proportion of Participants',
                     limits = c(0,1),
                     expand = c(0,0)) +
  scale_x_discrete('Subtype',
                   labels = c('A' = 'A',
                              'B' = 'B',
                              'C' = 'C', 
                              'D' = 'D',
                              'CRFO1_AE' = 'CRFO1_AE',
                              'CRFO2_AG' = 'CRFO2_AG',
                              'recombinant'= 'Recombinant',
                              'unknown' = 'Unknown') %>%
                     str_wrap(width = 17),
                   expand = c(0,0))+
  my_theme+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.y = element_blank())


# Fig 2E: Sampling Delay (inferred from feibig/seropositivty/timing of sampling relative to infection)
plt_2e <- ggplot(df, aes(x = sampling.delay_))+ 
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_)), 
               y = (..count..)/sum(..count..)))+
  scale_fill_brewer("Founder Multiplicity",
                    palette = 'YlGn',
                    labels = labs) +
  scale_y_continuous('Proportion of Participants',
                     limits = c(0,1),
                     expand = c(0,0)) +
  scale_x_discrete('Sampling Delay',
                   labels = c('<21' = '<21 Days', 
                              '>21' = '>21 Days',
                              'Unknown' = 'Unknown') %>% 
                     str_to_title(),
                   expand = c(0,0)) +
  my_theme+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.y = element_blank())


# Fig 2F: Number of consensus sequences analysed
numseqs_df <- GetNumSeqs(df)

plt_2f <- ggplot(numseqs_df , aes(sequencing.number_))+
  geom_histogram(aes(fill = forcats::fct_rev(factor(multiple.founders_)), 
                     y = (..count..)/sum(..count..)), 
                 binwidth = 3)+
  scale_fill_brewer("Founder Multiplicity",
                    palette = 'YlGn',
                    labels = labs) +
  scale_y_continuous('Proportion of Participants',
                     limits = c(0,1),
                     expand = c(0,0)) +
  scale_x_continuous('Number of Consensus Sequences',
                     limits = c(0,160),
                     expand = c(0,0))+
  my_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.y = element_blank())


# Fig 2G: Genomic region analysed
plt_2g <- ggplot(df, aes(x = sequencing.gene_))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_)),
               y = (..count..)/sum(..count..)))+
  scale_fill_brewer("Founder Multiplicity",
                    palette = 'YlGn',
                    labels = labs) +
  scale_y_continuous('Proportion of Participants',
                     limits = c(0,1),
                     expand = c(0,0)) +
  scale_x_discrete('Genomic Region Analysed',
                   labels = c('env' = 'Env',
                              'gag'= 'Gag',
                              'pol'='Pol',
                              'whole.genome' = 'NFLG'),
                   expand = c(0,0))+
  my_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.y = element_blank())


# Fig 2H =  Alignment Length
plt_2h <- ggplot(df, aes(alignment.length_))+
  geom_histogram(aes(fill = forcats::fct_rev(factor(multiple.founders_)), 
                     y = (..count..)/sum(..count..)), binwidth = 100)+
  scale_fill_brewer("Founder Multiplicity",
                    palette = 'YlGn',
                    labels = labs) +
  scale_y_continuous('Proportion of Participants',
                     limits = c(0,1),
                     expand = c(0,0)) +
  scale_x_continuous('Alignment Length',
                     limits = c(0,9500),
                     expand = c(0,0)) +
  my_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.y = element_blank())


# Fig 2I =  SGA
df_sga <- df %>% mutate(sga = ifelse(sample.amplification_ == 'SGA', 'SGA', 'Not SGA'))

plt_2i <- ggplot(df_sga, aes(x = sga))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders_)), 
               y = (..count..)/sum(..count..)))+
  scale_fill_brewer("Founder Multiplicity",
                    palette = 'YlGn',
                    labels = labs) +
  scale_y_continuous('Proportion of Participants',
                     limits = c(0,1),
                     expand = c(0,0)) +
  scale_x_discrete('Amplification', expand = c(0,0)) +
  my_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.y = element_blank())


panel_2 <- cowplot::plot_grid(plt_2b + theme(legend.position= c(0.85,0.8), axis.title.x = element_text(margin = margin(t = 5))),
                              plt_2c + theme(legend.position="none", axis.title.x = element_text(margin = margin(t = 12))),
                              plt_2d + theme(legend.position="none", axis.title.x = element_text(margin = margin(t = 12))), 
                              plt_2e + theme(legend.position="none", axis.title.x = element_text(margin = margin(t = 21))),
                              plt_2f + theme(legend.position="none", axis.title.x = element_text(margin = margin(t = 9))),
                              plt_2g + theme(legend.position="none"),
                              plt_2h + theme(legend.position="none", axis.title.x = element_text(margin = margin(t = 8))),
                              plt_2i + theme(legend.position="none", axis.title.x = element_text(margin = margin(t = 0))),
                            ncol = 2, 
                            nrow = 4,
                            align = "hv", 
                            axis = "bt",
                            labels = 'AUTO') 


####################################### Figure 3 - compare pooling models #######################################
pooled_models <-  read.csv(paste('results',sep = '/', 'pooling_estsa2sa3sa4sa6sa7.csv'))

resampled_models <- read.csv(paste('results',sep = '/', 'pooling_sa5.csv'))

models <- c('Two-Step Binomial-Normal',
            'One-Step Binomial GLMM')

og_models <- cbind("model" = pooled_models[1:2,1], pooled_models[1:2,3:8] %>% round(digits = 3)) %>%
  arrange(., model)
library(ggridges)

plt_3 <- ggplot(resampled_models) + 
  geom_density_ridges(aes(x = estimate,
                          y = analysis,
                      
                          fill=analysis), stat="binline",  binwidth = 0.0005,   scale = 1,  colour = NA,)+
  geom_point(inherit.aes = FALSE, data = og_models, 
             aes(x= estimate, y = model, color = model),position = position_nudge(y = -0.05), size = 3) +
  geom_errorbarh(inherit.aes = FALSE, data = og_models, 
                 aes(xmin = estimate.lb, xmax = estimate.ub, y = model, color = model), 
                 height = 0.15, position = position_nudge(y = -0.05)) + 

  scale_x_continuous(name = "Probability of Multiple Variants",
                     limits=c(0.2,0.4), 
                     expand = c(0, 0))+
  
  scale_y_discrete(name = "Model", 
                   labels = c('onestep_bi_rand'='One-step Binomial GLMM', 
                              'twostep_binorm' = 'Two-step Binomial-Normal'
                              )%>%
                     str_wrap(width = 17)) +
  
  scale_fill_brewer(palette = 'YlGn')+
  scale_colour_brewer(palette = 'YlGn')+
  my_theme


  

####################################### Figure 4 - meta-regression #######################################
pred <-  read.csv('./results/multimetareg_preds.csv')

level_order <- c("PWID","MTC:notiming",
                 "MTC:PostP",
                 "MTC:IntraP",
                 "MTC:PreP",
                 'MSM',
                 'HSX:nodirection',
                 'HSX:FTM',
                 'HSX:MTF')
pred$covariate_level <- factor(pred$covariate_level, levels = level_order) 
pred$grp <- c('HSX','HSX','HSX','MSM', 'MTC', 'MTC', 'MTC', 'MTC', 'PWID') 

plt_4 <- ggplot() +
  geom_linerange(aes(y = covariate_level, 
                     xmin= conf.low, 
                     xmax= conf.high), 
                 data = pred)+
  geom_point(aes(x= predicted, 
                 y = covariate_level,
                 size = tabs[[2]]),
             shape = 18,
             data = pred,
             colour = '#78c679') +
  scale_size(range = c(2.5,7.5)) +
  scale_x_continuous("Probability of Multiple Variants",
                     expand = c(0,0))+
  scale_y_discrete(labels = c("PWID" = 'PWID', 
                              "MTC:IntraP" = 'Mother-to-child: intrapartum', 
                              "MTC:notiming" = 'Mother-to-child: undisclosed',
                              "MTC:PostP" = 'Mother-to-child: post-partum',
                              "MTC:PreP" = 'Mother-to-child: pre-partum',
                              'MSM' = 'MSM', 
                              'HSX:nodirection' = 'Heterosexual: undisclosed',
                              'HSX:MTF' = 'Male-to-female',
                              'HSX:FTM' = 'Female-to-male'),
                   name = element_blank()) + 
  coord_cartesian(xlim = c(0,0.6)) +
  my_theme
  