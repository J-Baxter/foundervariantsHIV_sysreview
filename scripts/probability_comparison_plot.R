
pred <-  read.csv('./results/multimetareg_preds_env.csv') %>% select(-X)
patel <- data.frame(route = c('IPV', 'RPV', 'IAI', 'PWID', 'RAI', 'MTC'),
                    x = c(4,8,11,63,138,2260),
                    x.lb = c(1,6,4,41,102,1700),
                    x.ub = c(14,11,28,92,186,2900),
                    covariate_level = c('HSX:FTM', 'HSX:MTF', 'MSM', 'PWID', 'MSM', 'MTC:notiming'))

plotdata <- dplyr::left_join(patel, pred, by = "covariate_level") %>%
  mutate(x = x/10000) %>%
  mutate(x.lb = x.lb/10000) %>%
  mutate(x.ub = x.ub/10000) %>%
  filter(route != 'MTC')

                       
setEPS()
postscript(file = 'comparison_plot.eps', height = 8, width = 12)
ggplot(plotdata) +
  geom_point(aes(x = x, y = predicted, color = covariate_level))+
  geom_linerange(aes(xmin = x.lb, xmax = x.ub, y = predicted, color = covariate_level ))+
  geom_linerange(aes(ymin = conf.low, ymax = conf.high, x = x, color = covariate_level)) +
  scale_x_log10(name = 'Probability of Acquisition', expand = c(0,0)) + 
  scale_y_continuous(name = 'Probability of Multiple Founders', expand = c(0,0), limits = c(0,.6)) +
  scale_color_discrete(name = 'Exposure')+
  coord_flip() + 
  theme_classic()+annotation_logticks()+
  theme(legend.position = 'bottom',
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)
)
dev.off()
                  