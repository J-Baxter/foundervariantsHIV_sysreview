###################################################################################################
###################################################################################################
# Visualisation for pooling and sensitivity analyses
###################################################################################################
###################################################################################################
# Dependencies
library(ggplot2)
library(ggsci)
library(kableExtra)
library(metafor)
library(dplyr)
source('./scripts/generalpurpose_funcs.R')

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
influence_df <- read.csv("./results_nophylosplit/pooling_sa1.csv") %>% arrange(., model)

pooled_models <-  read.csv('./results_nophylosplit/pooling_estsa2sa3sa4.csv')

resampled_models <- read.csv('./results_nophylosplit/pooling_boot.csv')

models <- c('Two-Step Binomial-Normal',
            'One-Step Binomial GLMM')

og_models <- cbind("model" = pooled_models[1:2,1], pooled_models[1:2,3:8] %>% round(digits = 3)) %>% arrange(., model)


###################################################################################################
###################################################################################################
# Panel: Table of estimate and tau for pooling

og_models_formatted <- cbind.data.frame("Estimate" = paste0(og_models$estimate, ' ' ,
                                                               '[' , og_models$estimate.lb ,' - ',
                                                               og_models$estimate.ub, ']'),
                                        og_models[,5:7],
                                        row.names = rev(models) ) 
og_models_formatted[is.na(og_models_formatted )] <- "-"

# Not hapy with formatting
tbl <- kbl(og_models_formatted , 
          booktabs = T,
           col.names = c('Estimate', '$$\\hat{\\tau}^2$$', "$$\\text{Q}$$", "$$\\text{I}^2$$"),
           escape = FALSE,
           align = 'c', 
           ) %>%
  add_header_above(c(" " = 1, "Probability of Multiple Founders" = 1, "Heterogeneity" = 3), bold = T, line = T) %>%
 kable_classic(full_width = F, html_font = 'arial')


###################################################################################################
###################################################################################################
# Panel: Sensitivity analyses (exclusion criteria and resampling)

# Exclusion critera dot and whisker plot
senseplot <- ggplot(pooled_models,
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
        axis.text = element_text(size = 9.5,  family = "sans"),
        legend.text = element_text(size = 9.5,  family = "sans"),
        axis.title = element_text(size = 11,  family = "sans"),
        legend.background = element_blank()#,
        #plot.margin = unit(c(2,4,2,1), "lines")
  )


boot_plt <- ggplot(resampled_models) + 
  geom_histogram(aes(x = estimate,
                 color=analysis, 
                 fill=analysis),
                 binwidth = 0.0005) + 
  
  #geom_vline(aes(xintercept=mean(estimate)),color="#DC0000B2", linetype="dashed", size=1) +
  geom_vline(aes(xintercept= estimate,
            color= model),
             linetype="dashed", 
             size=0.75, 
            data = og_models) +
  
  theme_bw() +
  
  scale_y_continuous(name = "Frequency" ,
                     limits=c(0,75),
                     expand = c(0, 0)) +
  
  scale_x_continuous(name = "Probability of Multiple Founders",
                     limits=c(0,0.5), 
                     expand = c(0, 0))+

  geom_rect(data = og_models,
            aes(ymin=0, ymax=Inf,
            xmin=estimate.lb, xmax=estimate.ub,fill = model),alpha = .1) +
  scale_color_npg(name = c(analysis = 'Model'),
                  labels = c(
    onestep_bi_rand = "One-Step GLMM",
    twostep_binorm = "Two-Step Binomial Normal"))+
  scale_fill_npg(name = c(analysis = 'Model'),
                 labels = c(
    onestep_bi_rand = "One-Step GLMM",
    twostep_binorm = "Two-Step Binomial Normal"))+
  
  theme(axis.text = element_text(size = 9.5,  family = "sans"),
        legend.text = element_text(size = 9.5,  family = "sans"),
        axis.title = element_text(size = 11,  family = "sans"),
        legend.position = c(0.82,0.9),
        legend.background = element_blank())



jpeg("./results/pooling_sa.jpeg" ,width = 5000, height = 2500, res = 380 ,units = "px", pointsize = 12)
cowplot::plot_grid(senseplot, 
                   boot_plt , ncol = 2,  rel_widths  = c(1,1) ,labels = "AUTO", align = 'h', axis = 'b', greedy = F)
dev.off()

ggsave("./results/figures4.pdf", width = 16, height = 10, units= 'in')
cowplot::plot_grid(senseplot, 
                   boot_plt , ncol = 2,  rel_widths  = c(1,1) ,labels = "AUTO", align = 'h', axis = 'b', greedy = F)
dev.off()

jpeg("./results/boot_sensitivity.jpeg", width = 10, height = 16, units= 'in')
cowplot::plot_grid(boot_plt ,sa5_plt,  ncol = 2,  rel_widths  = c(1,1) ,labels = "AUTO", align = 'h', axis = 'b', greedy = F)
dev.off()
###################################################################################################
# Panel: Sensitivity analyses (Faceted Influence Plots)
influence_df$trial <- gsub("_" , " ", influence_df$trial)

plt <- ggplot(influence_df,aes(x = trial , y = estimate) ) +
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

jpeg(filename = './results/pooling_sa1.jpeg', res = 380, width=3600, height=5000 , units = 'px', pointsize = 10)

plt

dev.off()

ggsave("./results/figure_S5.pdf", width = 10, height = 10, units= 'in')
plt
dev.off()

###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################