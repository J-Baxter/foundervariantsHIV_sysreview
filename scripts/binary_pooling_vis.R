###################################################################################################
###################################################################################################
# Visualisation for pooling and sensitivity analyses
###################################################################################################
###################################################################################################
# RUN FROM HERE #
# Dependencies
source('./scripts/load_packages.R')
source('./scripts/generalpurpose_funcs.R')
require(ggplot2)
require(ggsci)
require(grDevices)

# Plot for pseudo-bootstrap replicates of participant selection
PltBoot <- function(data, intercept, ci.lb, ci.ub){
  plt <- ggplot(data) + 
    geom_histogram(aes(x = estimate),
                   color="black", 
                   fill="grey96",
                   binwidth = 0.0005) + 
    
    #geom_vline(aes(xintercept=mean(estimate)),color="#DC0000B2", linetype="dashed", size=1) +
    geom_vline(aes(xintercept= intercept),
               color="#3CBB7577", #00A087B2
               linetype="dashed", 
               size=1) +
    
    theme_classic() +
    
    scale_y_continuous(name = "Frequency" ,
                       limits=c(0,75),
                       expand = c(0, 0)) +
    
    scale_x_continuous(name = "Probability of Multiple Founders",
                       limits=c(0,0.35), 
                       expand = c(0.01, 0.01)) +
    
    annotate( "rect",
              ymin=0, ymax=Inf,
              xmin=ci.lb, xmax=ci.ub
              ,alpha = .1,
              fill = "#3CBB7577") +
    
    theme(axis.text = element_text(size = 10.5,  family = "sans"),
          legend.text = element_text(size = 10.5,  family = "sans"),
          axis.title = element_text(size = 13,  family = "sans"))
  
  return(plt)
}


###################################################################################################
###################################################################################################
# Import data
results_dir = './results'
influence_df <- read.csv(paste(results_dir,sep = '/', "pooling_sa1.csv")) %>% arrange(., model)
influence_rg <- read.csv(paste(results_dir,sep = '/', "pooling_sa8.csv")) %>% arrange(., model)

pooled_models <-  read.csv(paste(results_dir,sep = '/', 'pooling_estsa2sa3sa4sa6sa7.csv'))

resampled_models <- read.csv(paste(results_dir,sep = '/', 'pooling_sa5.csv'))

models <- c('Two-Step Binomial-Normal',
            'One-Step Binomial GLMM')

og_models <- cbind("model" = pooled_models[1:2,1], pooled_models[1:2,3:8] %>% round(digits = 3)) %>% arrange(., model)


###################################################################################################
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
  theme_bw() + 
  
  coord_flip() +
  
  geom_linerange(aes(ymin=estimate.lb, 
                     ymax= estimate.ub, 
                     color = analysis), 
                 position = position_dodge(0.5)) +
  
  scale_color_npg(name = 'Analysis', labels = c(
    no_small = "Studies with (n<10) omitted",
    no_zeros = "Studies with (p=0) omitted",
    original = "Full analysis",
    sga_only = 'Only SGA sequences')) + 
  
  theme(legend.position = c(0.8,0.86),
        axis.text = element_text(size = 9.5),
        legend.text = element_text(size = 9.5),
        axis.title = element_text(size = 11),
        legend.background = element_blank()#,
        #plot.margin = unit(c(2,4,2,1), "lines")
  )

# S2B
# Resampled Models
figureS2_b  <- ggplot(resampled_models) + 
  geom_histogram(aes(x = estimate,
                 color=analysis, 
                 fill=analysis),
                 binwidth = 0.0005) + 
  
  geom_vline(aes(xintercept= estimate,
            color= model),
             linetype="dashed", 
             size=0.75, 
            data = og_models) +
  
  theme_bw() +
  
  scale_y_continuous(name = "Frequency" ,
                     limits=c(0,120),
                     expand = c(0, 0)) +
  
  scale_x_continuous(name = "Probability of Multiple Founders",
                     limits=c(0.2,0.4), 
                     expand = c(0, 0))+
  

  geom_rect(data = og_models,
            aes(ymin=0, ymax=Inf,
            xmin=estimate.lb, xmax=estimate.ub,fill = model),alpha = .1) + #refigure as normal dist (from mean and calculate sd from CI's)
  scale_color_npg(name = c(analysis = 'Model'),
                  labels = c(
    onestep_bi_rand = "One-Step GLMM",
    twostep_binorm = "Two-Step Binomial Normal"))+
  scale_fill_npg(name = c(analysis = 'Model'),
                 labels = c(
    onestep_bi_rand = "One-Step GLMM",
    twostep_binorm = "Two-Step Binomial Normal"))+
  
  theme(axis.text = element_text(size = 9.5),
        legend.text = element_text(size = 9.5),
        axis.title = element_text(size = 11),
        legend.position = c(0.82,0.9),
        legend.background = element_blank())

# S3
# Effect of number of genomes: Dot and whisker
# All, no extreme, no small, no high
analysis_order <- c('original',
                 'no_extreme', 
                 'smallgenomes', 
                 'largegenomes',
                 'gold_standard') %>% rev()

plots3_data <- pooled_models[c(1,2,9:16),]
plots3_data$analysis <- factor(plots6_data$analysis, levels = analysis_order) 

figureS3 <- ggplot(plots6_data ,
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
  theme_bw() + 
  
  coord_flip() +
  
  guides(colour = guide_legend(reverse=T))+
  
  geom_linerange(aes(ymin=estimate.lb, 
                     ymax= estimate.ub, 
                     color = analysis), 
                 position = position_dodge(0.5)) +
  
  scale_colour_npg(name = 'Analysis', labels = c(
    original = "Full analysis",
    gold_standard = "Restricted to 'gold-standard' methodology",
    no_extreme = 'Restricted to 11-28 genomes/patient',
    smallgenomes = "Restriced to <11 genomes/patient",
    largegenomes = "Restricted to >28 genomes/patient")) + 
  
  theme(legend.position = c(0.8,0.86),
        axis.text = element_text(size = 9.5),
        legend.text = element_text(size = 9.5),
        axis.title = element_text(size = 11),
        legend.background = element_blank()#,
        #plot.margin = unit(c(2,4,2,1), "lines")
  )


figureS2 <- cowplot::plot_grid(figureS2_a, 
                               figureS2_b, 
                               ncol = 2,  rel_widths  = c(1,1) ,labels = "AUTO", align = 'h', axis = 'b', greedy = F)



# Save to file (ggsave rather than setEPS() to preseve transparencies)
ggsave(paste(figs_dir,sep = '/', "figureS5.eps"), device=cairo_ps, width = 16, height = 10, units= 'in')
Sys.sleep(0.5)
figureS2
dev.off()

ggsave(paste(figs_dir,sep = '/', "figureS6.eps"), device=cairo_ps, width = 8, height = 8, units= 'in')
Sys.sleep(0.5)
figureS3
dev.off()


###################################################################################################
# Panel: Sensitivity analyses (Faceted Influence Plots)
influence_df$trial <- gsub("_" , " ", influence_df$trial)

figureS4 <- ggplot(influence_df,aes(x = trial , y = estimate) ) +
  geom_point() + 
  
  scale_y_continuous(limits=c(0,0.4),
                     expand = c(0, 0),
                     name = "Probability of Multiple Founders") +
  
  scale_x_discrete(labels = influence_df$trial) +
  
  theme_bw() + 
  
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
  
  scale_fill_npg()+
  scale_color_npg()+
  
  theme(axis.line.y =element_blank(),
    axis.title.y =element_blank(),
    axis.ticks.y=element_blank(),
    strip.text.x = element_text(face = "bold" , colour = 'black' , size = 10.5),
    panel.spacing.x = unit(0.6 , 'cm'),
    strip.background = element_rect(fill = NA, colour= NA),
    legend.position = 'none'
  )

ggsave(paste(figs_dir,sep = '/', "figureS7.eps"), device=cairo_ps, width = 8, height = 10, units= 'in')
Sys.sleep(0.2)
figureS4
dev.off()


figureS5 <- ggplot(influence_rg,aes(x = trial , y = estimate) ) +
  geom_point() + 
  
  scale_y_continuous(limits=c(0,0.4),
                     expand = c(0, 0),
                     name = "Probability of Multiple Founders") +
  
  scale_x_discrete(labels = influence_rg$trial) +
  
  theme_bw() + 
  
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
  
  scale_fill_npg()+
  scale_color_npg()+
  
  theme(axis.line.y =element_blank(),
        axis.title.y =element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.x = element_text(face = "bold" , colour = 'black' , size = 10.5),
        panel.spacing.x = unit(0.6 , 'cm'),
        strip.background = element_rect(fill = NA, colour= NA),
        legend.position = 'none'
  )

ggsave(paste(figs_dir,sep = '/', "figureS8.eps"), device=cairo_ps, width = 8, height = 5, units= 'in')
Sys.sleep(0.2)
figureS5
dev.off()
###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################