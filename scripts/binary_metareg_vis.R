###################################################################################################
###################################################################################################
# Visualisation for meta-regression
###################################################################################################
###################################################################################################
# RUN FROM HERE #
# Dependencies
source('./scripts/load_packages.R')
source('./scripts/generalpurpose_funcs.R')


PlotMetaReg <- function(data, var){
  mycols_founder <- c('#002366','#DC143C', '#000000')
  #c("#D6604D", "#4393C3", '#000000')
  data.subset <- data[which(data$covariate %in% var),]
  data.subset$arrow.start <-  ifelse(exp(data.subset$ci.ub) > 8, exp(data.subset$est), 100)
  data.subset$arrow.end <-  ifelse(exp(data.subset$ci.ub) > 8, 8, 100)
  
  
  plt <- ggplot(data.subset) +
    geom_point(aes(x= exp(est), 
                   y = fct_reorder(level, order),
                   col =  ifelse(exp(est)>1 & exp(ci.lb)>1, "A", ifelse(exp(est)<1 & exp(ci.ub)<1, "B",  'C')),
                   size = tabs),
                   shape = 18) +
    theme_bw() + 
    geom_linerange(aes(y = level, 
                       xmin= exp(ci.lb), 
                       xmax= exp(ci.ub), 
                       col = ifelse(exp(est)>1 & exp(ci.lb)>1, "A", ifelse(exp(est)<1 & exp(ci.ub)<1, "B",  'C'))))+
    scale_x_continuous(
                       expand = c(0,0), 
                       name = "Odds Ratio",
                       #trans = 'log10'
                       )+
    scale_colour_manual(values = setNames(c("#E64B35FF", "#4DBBD5FF", '#000000'), c('A',"B","C"))) +
    scale_size(range = c(2.5,7.5)) + 
    geom_vline(xintercept = 1, linetype = 'dashed')+
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = 'none',
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title.y = element_blank()) +
    coord_cartesian(xlim = c(0,8))+
    #facet_grid(ref~., drop = T, scales = 'free', margins = F, space = 'free')++
    geom_segment(aes(x = arrow.start  , xend = arrow.end, y = level, yend = level,
                     col =  ifelse(exp(est)>1 & exp(ci.lb)>1, "A", ifelse(exp(est)<1 & exp(ci.ub)<1, "B",  'C'))),
                 arrow = arrow(length = unit(0.25, 'cm')))
  
  return(plt)
}

FPlot <- function(data,plotname){
  
  p2 <- ggplot(data = data) + 
    geom_point(aes(x= var, y = est)) + 
    theme_classic() + 
    geom_linerange( aes(x = var, ymin=fix.lb,ymax=fix.ub))+
    scale_y_continuous(limits = c(0,0.75) , 
                       expand = c(0,0), 
                       name = "Probability of Multiple Founders")+
    scale_x_discrete(name = plotname, 
                     guide = guide_axis(angle = 50))+
    #coord_flip()+
    theme(
      axis.text = element_text(size = 10,  family = "sans"),
      legend.text = element_text(size = 10,  family = "sans"),
      axis.title.y = element_text(size = 11,  family = "sans"),
      axis.title.x = element_text(size = 12,  family = "sans"),
      plot.margin = unit(c(2,4,2,1), "lines")
    )
  
  return(p2)
}

OR2Percent <- function(log_odds_ratio){
  
  out <- ifelse(exp(log_odds_ratio) < 1, 
                (1/ exp(log_odds_ratio) - 1)*-100, 
                (exp(log_odds_ratio) - 1)*100  )
  
  return(out)
}

###################################################################################################
# Set directory and import results
fe <- read.csv('./results/multimetareg_fe_env.csv')
int <- read.csv('./results/multimetareg_int_env.csv')
re <- read.csv('./results/multimetareg_re_env.csv')
emm <- read.csv('./results/multimetareg_emm.csv')
pred <-  read.csv('./results/multimetareg_preds_env.csv')

# Import data
if (!dir.exists('data')){
  Retrieve('data.zip')
}else{
  Sys.sleep(0.2)
}

df <- read.csv("./data/meta_analysis_data.csv",
               na.strings = "NA",
               stringsAsFactors = T) %>% 
  formatDF(.,filter = c('reported.exposure','grouped.subtype','sequencing.gene', 'sampling.delay')) %>%
  filter(reported.exposure_ != 'unknown.exposure') %>%
  droplevels()

baseline.covar <- c("reported.exposure_", "grouped.method_", "grouped.subtype_","sequencing.gene_", "sampling.delay_",'alignment.bin_')
baseline.level <- c("HSX:MTF", "haplotype", "B" , "env" , "<21", 'NFLG')

df <- SetBaseline(df, baseline.covar, baseline.level)

tabs <- list(table(df$grouped.method_)[sort(names(table(df$grouped.method_)))] %>% as.integer(),
             table(df$reported.exposure_)[sort(names(table(df$reported.exposure_)))] %>% as.integer(),
             table(df$sampling.delay_)[sort(names(table(df$sampling.delay_)))] %>% as.integer(),
             table(df$sequencing.gene_)[sort(names(table(df$sequencing.gene_)))] %>% as.integer()) %>% unlist() 

#fe <- read.csv('./results/multimetareg_fe.csv')
fe <- cbind.data.frame(covariate = gsub('_', '', baseline.covar)[-c(3,6)], level = baseline.level[-c(3,6)], est = 1000) %>%
  plyr::rbind.fill(.,fe) %>%
  arrange(., covariate,level) %>%
  cbind.data.frame(., tabs) %>%
  arrange(., covariate,est)

index <- fe %>% 
  dplyr::count(covariate) %>%
  pull(var = n) %>% 
  sapply(., seq, from = 1, by = 1) %>% 
  unlist()

fe$order <- index

fe[fe == 1000] <- 0

###################################################################################################
# Plot Fixed Effects and CIs from selected model (output to eps)
plt.list <- lapply(c('reported.exposure','grouped.method', 'sequencing.gene', 'sampling.delay'), PlotMetaReg, data = fe) #selected.fe 
plt_grid1 <- cowplot::plot_grid(plt.list[[1]] + scale_y_discrete(labels = str_wrap(c('HSX:FTM' = 'HSX: female-to-male',
                                                                            "MTC:PreP" = 'MTC: pre-partum',
                                                                            "MTC:PostP" = 'MTC: post-partum',
                                                                            "MTC:notiming" = 'MTC: undisclosed',
                                                                            "MTC:IntraP" = 'MTC: intrapartum',
                                                                            'MSM' = 'MSM', 
                                                                            'HSX:nodirection' = 'HSX: undisclosed',
                                                                            "PWID" = 'PWID', 
                                                                            'HSX:MTF' = 'HSX: male-to-female'
                                                                            ), width = 13)),
                                plt.list[[2]]+scale_y_discrete(labels = str_wrap(c(
                                                                          "model" = 'Model', 
                                                                          "phylogenetic:R" = 'Phylogenetic: recipient',
                                                                          "phylogenetic:S&R" = 'Phylogenetic: paired',
                                                                          "distance" = 'Distance',
                                                                          "molecular" = 'Molecular',
                                                                          'haplotype' = 'Haplotype'), width = 13)),
                                plt.list[[3]]+scale_y_discrete(labels = str_wrap(c(
                                                                          "pol" = 'Pol',
                                                                          "whole.genome" = 'NFLG',
                                                                          "gag" = ' Gag', 
                                                                          "env" = 'Env'), width = 13)),
                                plt.list[[4]]+scale_y_discrete(labels = str_wrap(c(
                                                                          ">21" = '>21 Days',
                                                                          "unknown" = 'Unknown',
                                                                          "<21" = '<21 Days'), width = 13)),
                                align = 'hv', axis = 'b' , labels = "AUTO", ncol = 2,rel_widths = c(1,1), label_size = 12, vjust = 1.3)

## Change to EPS out
setEPS()
postscript("./results/metareg_ORplot2.eps", width = 10, height = 10)
plt_grid1 
dev.off()


###################################################################################################
# Estimates of frequency, stratified by route of transmission
level_order <- c("PWID",
                 "MTC:IntraP",
                 "MTC:notiming",
                 "MTC:PostP",
                 "MTC:PreP",
                 'MSM',
                 'HSX:nodirection',
                 'HSX:FTM',
                 'HSX:MTF')



pred$covariate_level <- factor(pred$covariate_level, levels = level_order) 
pred$grp <- c('HSX','HSX','HSX','MSM', 'MTC', 'MTC', 'MTC', 'MTC', 'PWID') 

figure_3A <- ggplot() +
  geom_point(aes(x= predicted, 
                 y = covariate_level,
                 size = tabs),
             shape = 18, data = pred) +
  scale_size(range = c(2.5,7.5)) +
  theme_bw() + 
  geom_linerange(aes(y = covariate_level, 
                     xmin= conf.low, 
                     xmax= conf.high), data = pred)+
  scale_x_continuous(
    expand = c(0,0), 
    name = "Probability of Multiple Founders"
    #labels = scales::percent,
    #sec.axis = dup_axis(breaks = 0)
  )+
  
  #geom_segment(aes(x= 0.5, y = 0.5, xend = 0.5, yend = 1), inherit.aes = F) +

  coord_cartesian(xlim = c(0,0.6)) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.text.x.top = element_blank(),
    legend.position = 'none',
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()) +
  #geom_rect(aes(ymin = Inf,
               # ymax =  6.5,
               # xmin = -Inf, 
               # xmax = Inf),
           # fill = 'grey',
           # alpha = 0.2) +
  #geom_rect(aes(ymin = 5.5,
              #  ymax =  1.5,
              #  xmin = -Inf, 
               # xmax = Inf),
            #fill = 'grey',
           # alpha = 0.2)+
  scale_y_discrete(labels = c("PWID" = 'PWID', 
                              "MTC:IntraP" = 'Mother-to-child: intrapartum', 
                              "MTC:notiming" = 'Mother-to-child: undisclosed',
                              "MTC:PostP" = 'Mother-to-child: post-partum',
                              "MTC:PreP" = 'Mother-to-child: pre-partum',
                              'MSM' = 'MSM', 
                              'HSX:nodirection' = 'Heterosexual: undisclosed',
                              'HSX:MTF' = 'Heterosexual: male-to-female',
                              'HSX:FTM' = 'Heterosexual: female-to-male'),
                   name = element_blank())
  
figure_3B <-  cowplot::plot_grid(
                                 plt.list[[2]]+scale_y_discrete(labels = str_wrap(c(
                                   "model" = 'Model', 
                                   "phylogenetic:R" = 'Phylogenetic: recipient',
                                   "phylogenetic:S&R" = 'Phylogenetic: paired',
                                   "distance" = 'Distance',
                                   "molecular" = 'Molecular',
                                   'haplotype' = 'Haplotype'), width = 13)),
                                 plt.list[[3]]+scale_y_discrete(labels = str_wrap(c(
                                   "pol" = 'Pol',
                                   "whole.genome" = 'NFLG',
                                   "gag" = ' Gag', 
                                   "env" = 'Env'), width = 13)),
                                 plt.list[[4]]+scale_y_discrete(labels = str_wrap(c(
                                   ">21" = '>21 Days',
                                   "unknown" = 'Unknown',
                                   "<21" = '<21 Days'), width = 13)),
                                 align = 'hv', axis = 'b' , labels = c('B', 'C', 'D'), nrow = 3 ,rel_heights = c(1,1,1), label_size = 12)
  


figure_3 <- cowplot::plot_grid(figure_3A,  figure_3B ,  labels = c('A'), ncol= 2, rel_widths = c(1.25,1))
postscript("./results/test_grid.eps", width = 10, height = 16)
test_grid <- plot_grid(figure_3A, plt.list[[3]]+scale_y_discrete(labels = str_wrap(c(
  "pol" = 'Pol',
  "whole.genome" = 'NFLG',
  "gag" = ' Gag', 
  "env" = 'Env'), width = 13)))
#EPS
setEPS()
postscript("./results/figure_3.eps", width = 10, height = 16)
#test_grid
figure_3

dev.off()

###################################################################################################
###################################################################################################
# differently formatted plot for OR for transmission (Figure S11A)
fe.subset <- fe[which(fe$covariate %in% 'reported.exposure'),]
fe.subset$arrow.start <-  ifelse(exp(fe.subset$ci.ub) > 5, exp(fe.subset$est), 100)
fe.subset$arrow.end <-  ifelse(exp(fe.subset$ci.ub) > 5, 5, 100)
level_order <- c("PWID",
                 "MTC:IntraP",
                 "MTC:notiming",
                 "MTC:PostP",
                 "MTC:PreP",
                 'MSM',
                 'HSX:nodirection',
                 'HSX:FTM',
                 'HSX:MTF')
fe.subset$level <- factor(fe.subset$level, levels = level_order) 

Figure_S10A <- ggplot() +
  geom_point(aes(x= exp(est), 
                 y = level,
                 col =  ifelse(exp(est)>1 & exp(ci.lb)>1, "A", ifelse(exp(est)<1 & exp(ci.ub)<1, "B",  'C'))),
             shape = 4, 
             size = 5,
             data = fe.subset) +
  theme_bw() + 
  geom_linerange(aes(y = level, 
                     xmin= exp(ci.lb), 
                     xmax= exp(ci.ub), 
                     col = ifelse(exp(est)>1 & exp(ci.lb)>1, "A", ifelse(exp(est)<1 & exp(ci.ub)<1, "B",  'C'))),
                 data = fe.subset)+
  scale_x_continuous(
    expand = c(0,0), 
    name = "Odds Ratio",
    #trans = 'log10'
  )+
  scale_y_discrete(labels = c('HSX:MTF' = 'Heterosexual: male-to-female',
                              'HSX:FTM' = 'Heterosexual: female-to-male',
                              'HSX:nodirection' = 'HSX: undisclosed',
                              'MSM' = 'MSM', 
                              "MTC:PreP" = 'Mother-to-child: pre-partum',
                              "MTC:PostP" = 'Mother-to-child: post-partum',
                              "MTC:notiming" = 'Mother-to-child: undisclosed',
                              "MTC:IntraP" = 'Mother-to-child: intrapartum', 
                              "PWID" = 'PWID'), drop = FALSE)+
  scale_colour_manual(values = setNames(c("#E64B35FF", "#4DBBD5FF", '#000000'), c('A',"B","C"))) +
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title.y = element_blank()) +
  coord_cartesian(xlim = c(0,5))+
  #facet_grid(ref~., drop = T, scales = 'free', margins = F, space = 'free')++
  geom_segment(aes(x = arrow.start  , xend = arrow.end, y = level, yend = level,
                   col =  ifelse(exp(est)>1 & exp(ci.lb)>1, "A", ifelse(exp(est)<1 & exp(ci.ub)<1, "B",  'C'))),
               arrow = arrow(length = unit(0.5, 'cm')),
               data = fe.subset)

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
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'bottom',
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title.y = element_blank()) +
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
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'bottom',
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title.y = element_blank()) +
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
             size = 3, 
             position = position_dodge2(width = 0.7)) +
  theme_bw() + 
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
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0.8,0.93),
    legend.background = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.placement = 'outside',
    axis.title.y = element_blank(),
    strip.text.y = element_blank(),
    strip.background.x = element_blank()) +
  coord_cartesian(xlim = c(0,6.5))+
  facet_grid(covariate ~ .,  scales = 'free_y', space = 'free_y', drop = T, switch = 'y' ) 
 

# Out to file
Figure_S10_L <- cowplot::plot_grid(Figure_S10A,
                                   Figure_S10B,
                                   Figure_S10C,
                                   labels = 'AUTO', nrow = 3, align = 'v', axis = 'l', rel_heights = c(1,1) ,vjust = 1)

Figure_S10 <- cowplot::plot_grid(Figure_S10_L, 
                                 Figure_S10D,
                                 ncol = 2, align = 'h', axis = 't',labels = c('','D'))

setEPS()
postscript("./results/figureS11.eps", width = 12, height = 12)
Sys.sleep(0.2)
Figure_S9
dev.off()


###################################################################################################
# END # 
###################################################################################################
###################################################################################################