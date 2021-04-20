###################################################################################################
###################################################################################################
# Visualisation for meta-regression
###################################################################################################
###################################################################################################
# Dependencies
library(ggplot2)
library(ggsci)
library(kableExtra)
library(metafor)
library(dplyr)
library(magick)
PlotMetaReg <- function(data, var){
  mycols_founder <- c(RColorBrewer::brewer.pal(name = 'RdBu', n = 8)[c(2,7)], '#000000')
  data.subset <- data[which(data$covariate %in% var),]
  
  plt <- ggplot(data.subset) +
    geom_point(aes(x= exp(est), 
                   y = reorder(level,est) ,
                   col =  ifelse(exp(est)>1 & exp(ci.lb)>1, "A", ifelse(exp(est)<1 & exp(ci.ub)<1, "B",  'C'))),
                   shape = 4, 
                   size = 6) +
    theme_bw() + 
    geom_linerange(aes(y = level, 
                       xmin= exp(ci.lb), 
                       xmax= exp(ci.ub), 
                       col = ifelse(exp(est)>1 & exp(ci.lb)>1, "A", ifelse(exp(est)<1 & exp(ci.ub)<1, "B",  'C'))))+
    scale_x_continuous(#limits = c(0.0001,20),
                       expand = c(0,0), 
                       name = "Odds Ratio")+
    scale_colour_manual(values = setNames(c("#D6604D", "#4393C3", '#000000'), c('A',"B","C"))) +
    geom_vline(xintercept = 1, linetype = 'dashed')+
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = 'none',
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title.y = element_blank())
  
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

###################################################################################################
# Set directory and import results
setwd("../results")

fe <- read.csv('fixef_modelbuild_fe.csv')
int <- read.csv('fixef_modelbuild_int.csv')
re <- read.csv('fixef_modelbuild_re.csv')
model_selected.effectstruct <-"Reported Exposure + Grouped Method + Sequencing Gene + Sampling Delay" #GetName(model_selected.form, effects = 'fixed') #Requires var from binary_metareg.R

selected.fe <- fe[which(fe$analysis %in% model_selected.effectstruct),]

selected.int <- int[which(int$analysis %in% model_selected.effectstruct),]
selected.re <- re[which(re$analysis %in% model_selected.effectstruct),]


###################################################################################################
# Plot Fixed Effects and CIs from selected model (output to jpeg)
plt.list <- lapply(c('reported.exposure','grouped.method', 'sequencing.gene', 'sampling.delay'), PlotMetaReg, data = selected.fe) #selected.fe 
plt_grid1 <- cowplot::plot_grid(plt.list[[2]]+scale_y_discrete(labels = c("model" = 'Model', 
                                                                          "phylogenetic" = 'Phylogenetic',
                                                                          "distance" = 'Distance',
                                                                          "molecular" = 'Molecular')),
                                plt.list[[3]]+scale_y_discrete(labels = c("gag" = ' Gag', 
                                                                          "env" = 'Env',
                                                                          "pol" = 'Pol')),
                                plt.list[[4]]+scale_y_discrete(labels = c("unknown" = 'Unknown', 
                                                                          ">21" = '>21 Days')),
                                align = 'h', axis = 'b' , labels = c('B', 'C', 'D'), ncol = 3,rel_widths = c(1,1,1))

plt_grid2 <- cowplot::plot_grid(plt.list[[1]] + scale_y_discrete(labels = c("PWID" = 'PWID', 
                                                           "MTC:PreP" = 'MTC:Pre-Partum',
                                                           "MTC:PostP" = 'MTC:Post-Partum',
                                                           "MTC:notiming" = 'MTC:No Timing',
                                                           "MTC:IntraP" = 'MTC:Intra-Partum', 
                                                           'MSM' = 'MSM', 
                                                           'HSX:nodirection' = 'HSX:No Direction')),
                                plt_grid1 ,  labels = c('A'), nrow = 2, rel_heights = c(2,1))

jpeg(filename = 'metareg_plot.jpeg', width = 3000, height = 4000, res = 380 ,units = "px", pointsize = 12)

plt_grid2

dev.off()


###################################################################################################
# SA7 Plot
sa7_fe <- read.csv('s7_fe.csv', stringsAsFactors = F)
sa7_fe$delay.status <- base::strsplit(sa7_fe$analysis, '[_ & .]') %>%
  sapply(., "[[", 2)
sa7_fe$repeat.status <- base::strsplit(sa7_fe$analysis, '[.]') %>%
  sapply(., "[[", 2)

mycols_founder <- RColorBrewer::brewer.pal(name = 'RdBu', n = 8)[c(2,7)]
data.subset <- data[which(data$covariate %in% var),]

plt3.1 <- ggplot(sa7_fe[which(sa7_fe$delay.status %in% 'unknown'),]) +
  geom_point(aes(x= est, y = reorder(level,est) ,col =  est<0, shape = repeat.status), 
             size = 4,
             position = position_dodge(width = 1)) +
  theme_bw() + 
  geom_linerange(aes(y = level, xmin= ci.lb, xmax= ci.ub, col =  est<0, linetype = repeat.status),
                 position = position_dodge(width = 1)) +
  scale_x_continuous(limits = c(-3,3),
                     expand = c(0,0), 
                     name = "Log Odds Ratio")+
  scale_colour_manual(values = mycols_founder) +
  geom_vline(xintercept = 0, linetype = 'dashed')+
  facet_grid(covariate ~.,  scales = 'free_y', space = 'free_y', drop = T, switch = 'y') +
  labs(title = 'Unknown Delay') +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title.y = element_blank())


plt3.2 <- ggplot(sa7_fe[which(sa7_fe$delay.status %in% 'nounknown'),]) +
  geom_point(aes(x= est, y = reorder(level,est) ,col =  est<0, shape = repeat.status), 
             size = 4,
             position = position_dodge(width = 1)) +
  theme_bw() + 
  geom_linerange(aes(y = level, xmin= ci.lb, xmax= ci.ub, col =  est<0, linetype = repeat.status),
                 position = position_dodge(width = 1)) +
  scale_x_continuous(limits = c(-3,3),
                     expand = c(0,0), 
                     name = "Log Odds Ratio")+
  scale_colour_manual(values = mycols_founder) +
  geom_vline(xintercept = 0, linetype = 'dashed')+
  facet_grid(covariate ~.,  scales = 'free_y', space = 'free_y', drop = T, switch = 'y') +
  labs(title = 'No Unknowns') +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title.y = element_blank())

prow <- cowplot::plot_grid(plt3.1, plt3.2 ,ncol  =2)

jpeg(filename = 'sa7_plot.jpeg', width = 3000, height = 4000, res = 380 ,units = "px", pointsize = 12)

cowplot::plot_grid(prow, get_legend(plt3.2+guides(color = guide_legend(nrow = 1))+theme(legend.position = "bottom")), ncol = 1, rel_heights = c(1,0.1))

dev.off()

###################################################################################################
ggplot(funnel_data ) +
  geom_polygon(aes(x=x, y = y), data =  poldgpn ,fill = 'white', linetype = 'dashed' , color = 'black')  +
  geom_point( aes(y = se, x = b, colour = ), shape = 4, size = 3)+
  theme_classic() +
  scale_x_continuous(limits = c(-5 , 3), expand = c(0,0), name = 'Log Odds of Multiple Founders')+
  scale_y_reverse(limit=c(1.5,0),  expand = c(0,0), name = 'Standard Error') +
  
  geom_segment(aes(x=u, y =1.5, xend = u, yend=0)) +
  theme(panel.background = element_rect(fill = 'gray97' )) +
  scale_color_npg()


###################################################################################################
###################################################################################################
# Import data
setwd("./data")

###################################################################################################
# Plot and tabulate random effects structure selection

raneff_selection <- read.csv('raneff_selection.csv') 
raneff_ic <- raneff_selection[,c(1,2,3)] %>%
  reshape2::melt() %>%
  `colnames<-`(c('X' , 'criteria' , 'value'))

raneff_plot <- ggplot(raneff_selection) + 
  geom_point(aes(x = X, y = estimate))+
  geom_linerange(aes(x = X, ymin=estimate.lb, 
                     ymax= estimate.ub))+
  geom_line(aes(x = X,
                y = value/2500,
                color = criteria, 
                group = criteria), data = raneff_ic) +
  geom_point(aes(x = X,
                 y = value/2500,
                 color = criteria,
                 group = criteria), data = raneff_ic)+
  scale_y_continuous(name = 'Probability of Multiple Founders',
                     expand = c(0,0.02),
                     limits = c(0,1), 
                     sec.axis = sec_axis(~.*2500 , name = 'AIC/BIC'))+
  scale_x_discrete(name = 'Random Effects Structure',
                   labels = stringr::str_wrap(as.factor(raneff_selection$X), width = 30))+
  theme_classic() + 
  scale_color_npg()+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# Format dataframe for table
raneff_comb <- cbind.data.frame("Estimate" = paste0(round(raneff_selection$estimate, digits = 3), ' ' ,
                                                    '[' , round(raneff_selection$estimate.lb, digits = 3) ,' - ',
                                                    round(raneff_selection$estimate.ub, digits = 3), ']'),
                                row.names = raneff_selection$X )

raneff_tbl <- cbind.data.frame(raneff_comb, raneff_selection[,c(2,3)] %>% round(digits =2))

knitr::kable(raneff_tbl, 
                   booktabs = T,
                   col.names = c('Estimate', 'AIC', 'BIC'),
                   escape = FALSE,
                   align = 'c', 
                   linesep = c("\\addlinespace")) %>% kable_classic(full_width = F, html_font = 'arial')
#need to sort output to latex - current error with magick/ghostscript not seeing eye ot eye
###################################################################################################
# Plot and tabulate univariate fixed effects 
test_uni <- read.csv('fixef_univariate_fe.csv') 
plotnames <-test_uni$names
fixeff_uni.split <- split.data.frame(fixeff_uni , fixeff_uni$names)

subgroup_plotlist <- mapply(FPlot,subgroup_fe ,plotnames, SIMPLIFY = F )
subgroup_plot <- plot_grid(plotlist = subgroup_plotlist , labels = "AUTO" , align = 'hv', ncol = 2)

fixeff_uni_comb <- cbind.data.frame("Estimate" = paste0(transf.ilogit(test_uni$est) %>% round(digits = 3), ' ' ,
                                                    '[' , transf.ilogit(test_uni$ci.lb) %>% round(digits = 3) ,' - ',
                                                    transf.ilogit(test_uni$ci.ub) %>% round(digits = 3), ']'))

fixeff_uni_tbl <- cbind.data.frame(test_uni$analysis, test_uni$level, fixeff_uni_comb )
write.csv(fixeff_uni_tbl, file = 'fixef_univariate_formatted.csv')
knitr::kable(fixeff_uni_tbl, 
             booktabs = T,
             col.names = c('Estimate', 'AIC', 'BIC'),
             escape = FALSE,
             align = 'c', 
             linesep = c("\\addlinespace")) %>% kable_classic(full_width = F, html_font = 'arial')

fixedplot.df <- rbind.data.frame(selected.fe, selected.int)
###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################