# Thesis Plots for appendices of chapter 2

####################################### Dependencies #######################################
library(tidyverse)
library(ggplot2)
require(extrafont)
source('./scripts/generalpurpose_funcs.R')

#https://www.fontsquirrel.com/fonts/latin-modern-sans
font_import(pattern = "lmsans10*") 
loadfonts()

library(scales)
library(extrafont)


#https://www.fontsquirrel.com/fonts/latin-modern-sans
font_import(pattern = "lmsans10*") 
loadfonts()

my_theme <- theme_classic(base_family = "LM Sans 10")+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 8),
    axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
    strip.text  = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
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

results_dir = './results'
influence_df <- read.csv(paste(results_dir,sep = '/', "pooling_sa1.csv")) %>% arrange(., model)
influence_rg <- read.csv(paste(results_dir,sep = '/', "pooling_sa8.csv")) %>% arrange(., model)

pooled_models <-  read.csv(paste(results_dir,sep = '/', 'pooling_estsa2sa3sa4sa6sa7.csv'))

resampled_models <- read.csv(paste(results_dir,sep = '/', 'pooling_sa5.csv'))

models <- c('Two-Step Binomial-Normal',
            'One-Step Binomial GLMM')

og_models <- cbind("model" = pooled_models[1:2,1], pooled_models[1:2,3:8] %>% round(digits = 3)) %>% arrange(., model)




###################################################################################################

# Year of Publication ~ Frequency of Individuals, stacked by risk group
figureS1_a <- ggplot(df, aes(year_))+
  geom_histogram(aes(fill = riskgroup_, 
                     colour = riskgroup_), 
                 bins = 8) +
  scale_color_brewer(palette = 'YlGn') +
  scale_fill_brewer(palette = 'YlGn') +
  scale_y_continuous(limits = c(0,500), 
                     expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  my_theme + 
  theme(legend.position = 'bottom',
        axis.text.x=element_text(angle=45, hjust=1)) +
  
  labs(fill = 'Risk Group', 
       x = 'Year of Publication', 
       y = 'Participants')+
  guides(color = "none")+
  theme(legend.position = 'right')

# Year of Publication ~ Frequency of Individuals, stacked by method of enumeration
figureS1_b <- ggplot(df, aes(year_))+
  geom_histogram(aes(fill = grouped.method_, 
                     colour = grouped.method_), 
                 bins = 8) +
 # scale_color_manual(values = mycols_method[c(2,4,6,8,10,12)]) +
  scale_colour_brewer(palette = 'YlGn') + 
  scale_fill_brewer(#values = mycols_method[c(2,4,6,8,10,12)],
    palette = 'YlGn',
                    labels = c('Haplotype', 
                               'Distance', 
                               'Model',
                               'Molecular',
                               'Phylogenetic: Recipient Only',
                               'Phylogenetic: Source & Recipient')) +
  scale_y_continuous(limits = c(0,500), 
                     expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  
  labs(fill = 'Method', 
       x = 'Year of Publication', 
       y = 'Participants')+
  guides(color = "none")+
  my_theme + 
  theme(legend.position = 'right',
        axis.text.x=element_text(angle=45, hjust=1)) 

# Year of Publication ~ Frequency of Individuals, stacked by sequencing technology
df$sequencing.method_<- factor(df$sequencing.method_, levels = c('sanger_SGA', 
                                                                 '2G:illumina_miseq', 
                                                                 '2G:roche_454',
                                                                 '3G:PacBio_hifi',
                                                                 'sanger', 
                                                                 'sanger_precSGA',
                                                                 'unknown'))


figureS1_c <- ggplot(df, aes(year_))+
  geom_histogram(aes(fill = sequencing.method_, 
                     colour = sequencing.method_), 
                 bins = 8) +
  scale_color_brewer(palette = 'YlGn') +
  scale_fill_brewer(palette = 'YlGn',
                    labels = c('Sanger with SGA', 
                               'Illumina MiSeq', 
                               'Roche 454',
                               'PacBio HiFi',
                               'Sanger Only', 
                               'Sanger with SGA Precursor',
                               'Unknown'), na.value = 'lightgrey') +
  scale_y_continuous(limits = c(0,500), 
                     expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  
  labs(fill = 'Sequencing Technology', 
       x = 'Year of Publication', 
       y = 'Participants')+
  guides(color = "none") +# #fill = guide_legend(nrow = 3, byrow= TRUE))+
  my_theme + 
  theme(legend.position = 'right',
        axis.text.x=element_text(angle=45, hjust=1)) 

figureS1 <- cowplot::plot_grid(figureS1_a, 
                               figureS1_b , 
                               figureS1_c, 
                               ncol = 1, labels = 'AUTO', align = 'hv')

###################################################################################################

# Panel: Sensitivity analyses (exclusion criteria and resampling)

# S2A
# Studies with (n<10) omitted, Studies with (p=0) omitted", Full analysis",'Only SGA sequences', 
# 'Gold Standard Only'

figureS2_a <- ggplot(pooled_models[1:8,],
                     aes(x= forcats::fct_rev(model), y = estimate, color = analysis)) +
  
  geom_point( shape = 4, 
              size = 4,
              position = position_dodge(0.5)) + 
  
  scale_y_continuous(name = "Probability of Multiple Founders",
                     limits=c(0,.5),
                     expand = c(0.01, 0.01)) +
  
  scale_x_discrete(name = "Model", 
                   labels = c(
                     onestep_bi_rand = "GLMM",
                     twostep_binorm = "B-N"
                   )) +
  
  coord_flip() +
  
  geom_linerange(aes(ymin=estimate.lb, 
                     ymax= estimate.ub, 
                     color = analysis), 
                 position = position_dodge(0.5)) +
  
  scale_colour_brewer(palette = 'YlGn',
                      name = 'Analysis', labels = c(
                        no_small = "Studies with (n<10) omitted",
                        no_zeros = "Studies with (p=0) omitted",
                        original = "Full analysis",
                        sga_only = 'Only SGA sequences')) + 
  my_theme +
  
  theme(legend.position = c(0.8,0.86)
  )


# S3
# Effect of number of genomes: Dot and whisker
# All, no extreme, no small, no high
analysis_order <- c('original',
                    'no_extreme', 
                    'smallgenomes', 
                    'largegenomes',
                    'gold_standard') %>% rev()

plots3_data <- pooled_models[c(1,2,9:16),]
plots3_data$analysis <- factor(plots3_data$analysis, levels = analysis_order) 

figureS2_b <- ggplot(plots3_data ,
                   aes(x= forcats::fct_rev(model), y = estimate, color = analysis)) +
  
  geom_point( shape = 4, 
              size = 4,
              position = position_dodge(0.5)) + 
  
  scale_y_continuous(name = "Probability of Multiple Founders",
                     #limits=c(0,.5),
                     expand = c(0.01, 0.01)) +
  coord_cartesian(ylim = c(0,.5))+
  
  scale_x_discrete(name = "Model", 
                   labels = c(
                     onestep_bi_rand = "GLMM",
                     twostep_binorm = "B-N"
                   )) +
  coord_flip() +
  
  guides(colour = guide_legend(reverse=T))+
  
  geom_linerange(aes(ymin=estimate.lb, 
                     ymax= estimate.ub, 
                     color = analysis), 
                 position = position_dodge(0.5)) +
  
  scale_colour_brewer(name = 'Analysis', 
                      palette = 'YlGn',
                      labels = c(
    original = "Full analysis",
    gold_standard = "Restricted to 'gold-standard' methodology",
    no_extreme = 'Restricted to 11-28 genomes/patient',
    smallgenomes = "Restriced to <11 genomes/patient",
    largegenomes = "Restricted to >28 genomes/patient")) + 
  
  my_theme +
  
  theme(legend.position = c(0.8,0.86)
  )


figure_s2 <- cowplot::plot_grid(figureS2_a, figureS2_b, align = 'hv', ncol = 1, nrow = 2, labels = 'AUTO')


###################################################################################################
# Panel: Sensitivity analyses (Faceted Influence Plots)
influence_df$trial <- gsub("_" , " ", influence_df$trial)

figureS4 <- ggplot(influence_df,aes(x = trial , y = estimate) ) +
  geom_point() + 
  
  scale_y_continuous(limits=c(0,0.4),
                     expand = c(0, 0),
                     name = "Probability of Multiple Founders") +
  
  scale_x_discrete(labels = influence_df$trial) +
  
  coord_flip() +
  
  facet_grid(cols = vars(model),
             labeller = labeller(model = c(
               onestep_bi_rand = "One-Step GLMM",
               twostep_binorm = "Two-Step Binomial Normal"))) +
  
  geom_errorbar(data =influence_df, 
                aes(x = trial,
                    ymin=ci.lb,
                    ymax=ci.ub)) +
  
  geom_hline(data = og_models, 
             aes(yintercept = estimate, color = model),
             linetype="dashed", 
             
             size=0.75) +
  
  geom_rect(data = og_models, 
            inherit.aes = FALSE,
            aes(ymin=estimate.lb, ymax=estimate.ub, 
                xmin = -Inf, xmax = Inf,fill = model),
            alpha = 0.1) +
  
  scale_fill_brewer(palette = 'YlGn')+
  scale_color_brewer(palette = 'YlGn')+
  
  theme_bw(base_family = "LM Sans 10")+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    strip.background = element_blank(),
    axis.title.y =element_blank(),
    axis.ticks.y=element_blank(),
    strip.text.x = element_text(face = "bold" , colour = 'black' , size = 8),
    panel.spacing.x = unit(0.6 , 'cm'),
    legend.position = 'none'
  )




figureS5 <- ggplot(influence_rg,aes(x = trial , y = estimate) ) +
  geom_point() + 
  
  scale_y_continuous(limits=c(0,0.4),
                     expand = c(0, 0),
                     name = "Probability of Multiple Founders") +
  
  scale_x_discrete(labels = influence_rg$trial) +
  
  coord_flip() +
  
  facet_grid(cols = vars(model),
             labeller = labeller(model = c(
               onestep_bi_rand = "One-Step GLMM",
               twostep_binorm = "Two-Step Binomial Normal"))) +
  
  geom_errorbar(data =influence_rg, 
                aes(x = trial,
                    ymin=ci.lb,
                    ymax=ci.ub),
                width = 0.2) +
  
  geom_hline(data = og_models, 
             aes(yintercept = estimate, color = model),
             linetype="dashed", 
             
             size=0.75) +
  
  geom_rect(data = og_models, 
            inherit.aes = FALSE,
            aes(ymin=estimate.lb, ymax=estimate.ub, 
                xmin = -Inf, xmax = Inf,fill = model),
            alpha = 0.1) +
  
  scale_fill_brewer(palette = 'YlGn')+
  scale_color_brewer(palette = 'YlGn')+
  
  theme_bw(base_family = "LM Sans 10")+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    strip.background = element_blank(),
    axis.title.y =element_blank(),
    axis.ticks.y=element_blank(),
    strip.text.x = element_text(face = "bold" , colour = 'black' , size = 8),
    panel.spacing.x = unit(0.6 , 'cm'),
    legend.position = 'none'
  )





# Save to file (ggsave rather than setEPS() to preseve transparencies)
ggsave('figureS1.eps', device=cairo_ps,  height = 210, width = 150, units = 'mm')
Sys.sleep(0.5)
figureS1
dev.off()


ggsave('figureS2.eps', device=cairo_ps,  height = 150, width = 150, units = 'mm')
Sys.sleep(0.5)
figureS2_a
dev.off()

Sys.sleep(0.5)

ggsave('figureS3.eps', device=cairo_ps,  height = 150, width = 150, units = 'mm')
Sys.sleep(0.5)
figureS2_b
dev.off()

Sys.sleep(0.5)

ggsave('figureS4.eps', device=cairo_ps,  height = 220, width = 150, units = 'mm')
Sys.sleep(0.5)
figureS4
dev.off()

Sys.sleep(0.5)

ggsave('figureS5.eps', device=cairo_ps,  height = 110, width = 150, units = 'mm')
Sys.sleep(0.5)
figureS5
dev.off()



###################################################################################################
# Sensitivity plots (For brevity, only the transmission covariates are plotted)
# S11B - LOOCV
sa1 <- read.csv('./results/multimetareg_s1.csv', stringsAsFactors = F) %>%
  filter(grepl('reported.exposure',X))

level_order <- c("reported.exposure_PWID",
                 "reported.exposure_MTC:IntraP",
                 "reported.exposure_MTC:notiming",
                 "reported.exposure_MTC:PostP",
                 "reported.exposure_MTC:PreP",
                 'reported.exposure_MSM',
                 'reported.exposure_HSX:nodirection',
                 'reported.exposure_HSX:FTM',
                 'reported.exposure_HSX:MTF')

sa1$level <- factor(gsub('[[:digit:]]' , '' , sa1$X), levels = level_order) 
colnames(sa1) <- c('X', 'trial', 'est', 'se', 'z.val', 'p.val', 'level')

Figure_S10B <-  ggplot() +
  geom_point(aes(x = exp(est), y =  level ), position = position_jitter(), data = sa1)+
  geom_point(aes(x = 1, y = 9))+
  theme_bw() + 
  scale_x_continuous(
    expand = c(0,0), 
    name = "Odds Ratio"
    #trans = 'log10'
  )+
  scale_y_discrete(labels = c('reported.exposure_HSX:MTF' = 'Heterosexual: male-to-female',
                              'reported.exposure_HSX:FTM' = 'Heterosexual: female-to-male',
                              'reported.exposure_HSX:nodirection' = 'HSX: undisclosed',
                              'reported.exposure_MSM' = 'MSM', 
                              "reported.exposure_MTC:PreP" = 'Mother-to-child: pre-partum',
                              "reported.exposure_MTC:PostP" = 'Mother-to-child: post-partum',
                              "reported.exposure_MTC:notiming" = 'Mother-to-child: undisclosed',
                              "reported.exposure_MTC:IntraP" = 'Mother-to-child: intrapartum', 
                              "reported.exposure_PWID" = 'PWID'), drop = FALSE)+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_family = "LM Sans 10")+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text  = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.position = 'bottom',
    legend.background = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank(),
    strip.placement = 'outside',
    axis.title.y = element_blank(),
    strip.text.y = element_blank() ) +
  coord_cartesian(xlim = c(0,5))


###################################################################################################
# S11C - boot resample
sa5_plotdata <- read.csv('./results/multimetareg_s5.csv', stringsAsFactors = F)%>%
  filter(grepl('reported.exposure',X))

level_order <- c("reported.exposure_PWID",
                 "reported.exposure_MTC:IntraP",
                 "reported.exposure_MTC:notiming",
                 "reported.exposure_MTC:PostP",
                 "reported.exposure_MTC:PreP",
                 'reported.exposure_MSM',
                 'reported.exposure_HSX:nodirection',
                 'reported.exposure_HSX:FTM',
                 'reported.exposure_HSX:MTF')

sa5_plotdata$level <- factor(gsub('[[:digit:]]' , '' , sa5_plotdata$X), levels = level_order) 
colnames(sa5_plotdata) <- c('X',  'est', 'se', 'z.val', 'p.val', 'rep','level')

Figure_S10C <- ggplot() +
  geom_point(aes(x = exp(est), y =  level ), position = position_jitter(), data = sa5_plotdata)+
  geom_point(aes(x = 1, y = 9))+
  scale_x_continuous(
    expand = c(0,0), 
    name = "Odds Ratio"
    #trans = 'log10'
  )+
  scale_y_discrete(labels = c('reported.exposure_HSX:MTF' = 'Heterosexual: male-to-female',
                              'reported.exposure_HSX:FTM' = 'Heterosexual: female-to-male',
                              'reported.exposure_HSX:nodirection' = 'HSX: undisclosed',
                              'reported.exposure_MSM' = 'MSM', 
                              "reported.exposure_MTC:PreP" = 'Mother-to-child: pre-partum',
                              "reported.exposure_MTC:PostP" = 'Mother-to-child: post-partum',
                              "reported.exposure_MTC:notiming" = 'Mother-to-child: undisclosed',
                              "reported.exposure_MTC:IntraP" = 'Mother-to-child: intrapartum', 
                              "reported.exposure_PWID" = 'PWID'), drop = FALSE)+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_family = "LM Sans 10")+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text  = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.position = 'bottom',
    legend.background = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank(),
    strip.placement = 'outside',
    axis.title.y = element_blank(),
    strip.text.y = element_blank() ) +
  coord_cartesian(xlim = c(0,5))



###################################################################################################
# S2-4 # to revisit
sa234_fe <- read.csv('./results/multimetareg_s2-4_fe.csv', stringsAsFactors = F) %>%filter(grepl('reported.exposure',covariate))
level_order <- c("PWID",
                 "MTC:IntraP",
                 "MTC:notiming",
                 "MTC:PostP",
                 "MTC:PreP",
                 'MSM',
                 'HSX:nodirection',
                 'HSX:FTM',
                 'HSX:MTF')
sa234_fe$level <- factor(sa234_fe$level, levels = level_order) 

Figure_S10D <- ggplot(sa234_fe) +
  geom_point(aes(x= exp(est), 
                 y = level,
                 col =  ifelse(exp(est)>1 & exp(ci.lb)>1, "A", ifelse(exp(est)<1 & exp(ci.ub)<1, "B",  'C')),
                 shape = analysis),
             position = position_dodge2(width = 0.7)) +
  scale_shape(name = 'Analysis', labels = c('No Small', 'No Zero', 'SGA Only'))+
  geom_linerange(aes(y = level, 
                     xmin= exp(ci.lb), 
                     xmax= exp(ci.ub), 
                     col = ifelse(exp(est)>1 & exp(ci.lb)>1, "A", ifelse(exp(est)<1 & exp(ci.ub)<1, "B",  'C'))),
                 position = position_dodge2(width = 0.7))+
  scale_x_continuous(
    expand = c(0,0), 
    name = "Odds Ratio"
    #trans = 'log10'
  )+
  scale_y_discrete(labels = c("PWID" = 'PWID', 
                              "MTC:PreP" = 'Mother-to-child: pre-partum',
                              "MTC:PostP" = 'Mother-to-child: post-partum',
                              "MTC:notiming" = 'Mother-to-child: undisclosed',
                              "MTC:IntraP" = 'Mother-to-child: intrapartum', 
                              'MSM' = 'MSM', 
                              'HSX:nodirection' = 'HSX: undisclosed',
                              'HSX:MTF' = 'Heterosexual: male-to-female',
                              'HSX:FTM' = 'Heterosexual: female-to-male',
                              'haplotype' = 'Haplotype', 
                              "model" = 'Model', 
                              "phylogenetic" = 'Phylogenetic',
                              "distance" = 'Distance',
                              "molecular" = 'Molecular',
                              "whole.genome" = 'NFLG',
                              "gag" = ' Gag', 
                              "env" = 'Env',
                              "pol" = 'Pol',
                              "unknown" = 'Unknown Delay', 
                              ">21" = '>21 Days',
                              "<21" = '<21 Days'))+
  scale_colour_manual(values = setNames(c("#E64B35FF", "#4DBBD5FF", '#000000'), c('A',"B","C")),guide = NULL) +
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_bw(base_family = "LM Sans 10")+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text  = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.position = c(0.8,0.93),
    legend.background = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank(),
    strip.placement = 'outside',
    axis.title.y = element_blank(),
    strip.text.y = element_blank() ) + 
  coord_cartesian(xlim = c(0,6.5)) +
  facet_grid(covariate ~ .,  scales = 'free_y', space = 'free_y', drop = T, switch = 'y' ) 


# Out to file
Figure_S10  <- cowplot::plot_grid(Figure_S10B,
                                   Figure_S10C,
                                   labels = 'AUTO', nrow = 3, align = 'v', axis = 'l', rel_heights = c(1,1) ,vjust = 1)


ggsave('figureS10.eps', device=cairo_ps,  height = 240, width = 150, units = 'mm')
Sys.sleep(0.5)
Figure_S10
dev.off()


ggsave('figureS11.eps', device=cairo_ps,  height = 200, width = 150, units = 'mm')
Sys.sleep(0.5)
Figure_S10D
dev.off()


###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################