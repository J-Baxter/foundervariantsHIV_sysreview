# Simulating the effect of sampling density on the observed number of founder lineages

# Dependencies
library(actuar)
library(ggplot2)
library(cowplot)
library(tidyverse)


# Dual founder infection with variable ratio between major/minor variant
Binom <-function(sample, prop_major){
  prop_minor <- 100 - prop_major 
  P_same <- (choose(prop_major,sample)+choose(prop_minor,sample))/choose(100,sample)
  P <- 1-P_same
  return(P)
}


# Zero-truncated binomial with varying probability (p) of selecting any given variant
ZT_Binom <- function(sample, p = 0.05){
  require(actuar)
  
  prob_1 <- dztbinom(1, sample, prob = p, log = FALSE)
  out <- 1-  prob_1
  return(out)
}


set.seed(4472)
sample_size <- 1:30


# Binomial 
major_props <-  seq(50,95, by = 5)
binom_probs <- sapply(sample_size, function(x) sapply(maj_props, Binom, sample = x)) %>% 
  as.data.frame() %>% 
  gather(key = 'sample_size', value = 'prob_both') %>% 
  cbind.data.frame(maj_proportion = seq(95,50, by = -5)) %>% 
  mutate(sample_size = replace(sample_size, values = rep(1:30, each = 10)))

binom_probs$sample_size <- as.numeric(binom_probs$sample_size)

plt_a <- ggplot(data = binom_probs,
                aes(x = sample_size, y =  prob_both, group= maj_proportion , color = maj_proportion)) + 
  geom_line() +
  scale_x_continuous(name = 'Sample Size', expand = c(0,0))+
  scale_y_continuous(name = 'P[Both Lineages Sampled]', expand = c(0,0)) +
  scale_color_viridis_c(name = 'Proportion of Major Variant (%)', direction = -1) +
  theme_classic() +
  theme(legend.position = 'none')


# Zero-truncated Binomial
p <- seq(0.05,0.8,by = 0.05)

zt_probs <- sapply(sample_size, function(x) sapply(p, ZT_Binom, sample = x)) %>% 
  as.data.frame() %>% 
  gather(key = 'sample_size', value = 'prob_more') %>% 
  cbind.data.frame(p = seq(0.8,0.05,by = -0.05)) %>% 
  mutate(sample_size = replace(sample_size, values = rep(1:30, each = 16)))

zt_probs$sample_size <- as.numeric(zt_probs$sample_size)

plt_b <- ggplot(data =zt_probs,
                aes(x = sample_size, y =  prob_more, group= p , color = p)) + 
  geom_line() +
  scale_x_continuous(name = 'Sample Size', expand = c(0,0))+
  scale_y_continuous(name = 'P[Finding >1 Lineages]', expand = c(0,0)) +
  scale_color_viridis_c(name = 'P[X=x]', direction = -1) +
  theme_classic()+
  theme(legend.position = 'none')


# Format output
fig <- plot_grid(plt_a + theme(legend.position = c(0.8, 0.25)), 
                 plt_b + theme(legend.position = c(0.8, 0.25)), 
                 align = 'hv', labels = 'AUTO')


# Save to file
setEPS()
postscript("./results/samplingprobs.eps", width = 16, height = 10)
fig
dev.off()

## END ##