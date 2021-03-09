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
fixeff_uni <- read.csv('fixeff_uni.csv') 
plotnames <-fixeff_uni$names
fixeff_uni.split <- split.data.frame(fixeff_uni , fixeff_uni$names)

subgroup_plotlist <- mapply(FPlot,subgroup_fe ,plotnames, SIMPLIFY = F )
subgroup_plot <- plot_grid(plotlist = subgroup_plotlist , labels = "AUTO" , align = 'hv', ncol = 2)

fixeff_uni_comb <- cbind.data.frame("Estimate" = paste0(round(fixeff_uni$est, digits = 3), ' ' ,
                                                    '[' , round(fixeff_uni$fix.lb, digits = 3) ,' - ',
                                                    round(fixeff_uni$fix.ub, digits = 3), ']'))

fixeff_uni_tbl <- cbind.data.frame(fixeff_uni$names, fixeff_uni$var, fixeff_uni_comb)

knitr::kable(raneff_tbl_df, 
             booktabs = T,
             col.names = c('Estimate', 'AIC', 'BIC'),
             escape = FALSE,
             align = 'c', 
             linesep = c("\\addlinespace")) %>% kable_classic(full_width = F, html_font = 'arial')

###################################################################################################
###################################################################################################
# END # 
###################################################################################################
###################################################################################################