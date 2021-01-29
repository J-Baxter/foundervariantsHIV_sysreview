###################################################################################################
###################################################################################################
#Figures describing covariate distirbution accompanying Inferring the multiplicity of founder variants
# initiating HIV-1 infection: a systematic review and IPD meta-analysis
# 1. Panel of bar plots (and one histogram) plotting frequency against covariate, stacking bars by 
# binary classification of multiplicity
# 2. 
###################################################################################################
###################################################################################################

# Dependencies
library(ggplot2)
library(ggsci)
library(dplyr)
library(forcats)
library(tidyr)
library(RColorBrewer)
library(cowplot)
source('generalpurpose_funcs.R')


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
GetNumSeqs <- function(df){
  # drop NA and NGS first 
  numseqs <- df[!df$sequencing.number %in% c('NGS', NA),]
  numseqs$sequencing.number <- as.numeric(numseqs_df$sequencing.number)
  
  return(numseqs)
}


# Group transmission variable at higher level (remove direction)
split_transission <- function(x , catnames){
  t1 <- as.character(x[,1]) %>% 
    strsplit(. , '[:]') %>% do.call(rbind.data.frame, .) %>% 
    cbind.data.frame(., x[,2])
  
  t2 <- t1[,c(1,3)]
  freq <- t2 %>%
    group_by(t2[,1] , t2[,2]) %>%
    summarise(frequency = n())
  colnames(freq) <- catnames
  return(t2)
}


###################################################################################################
###################################################################################################
# Import data
setwd("./data")
df <- read.csv("data_master_11121.csv", na.strings = "NA") %>% formatDF()


# Set colour palettes 
mycols_founder <- RColorBrewer::brewer.pal(name = 'RdBu', n = 8)[c(2,7)]
nb.cols <- 12
mycols_method <- colorRampPalette(brewer.pal(10, "RdBu"))(nb.cols)


###################################################################################################
###################################################################################################
# Basic plots for panel 1 (four bar plots to panel displaying frequencies of grouped method, 
# seropositivity, subtype and number of sequences)

# 1.1 Method
p1.method <- ggplot(df, aes(grouped.method))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders)), colour = forcats::fct_rev(factor(multiple.founders))))+
  scale_color_manual(values = mycols_founder)+
  scale_fill_manual(values = mycols_founder)+
  ylim(0,1000)+
  theme_classic()+
  xlab('Grouped Method')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")

# 1.2 Seroconversion (poss stack infant and NA together)
p1.seropos <- ggplot(df, aes(participant.seropositivity))+ 
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders)), colour = forcats::fct_rev(factor(multiple.founders))))+
  scale_color_manual(values = mycols_founder)+
  scale_fill_manual(values = mycols_founder)+
  ylim(0,1000)+
  theme_classic()+
  xlab('Seropositivity')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")

# 1.3 Subtype
p1.subtype <- ggplot(df, aes(grouped.subtype))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders)), colour = forcats::fct_rev(factor(multiple.founders))))+
  scale_color_manual(values = mycols_founder)+
  scale_fill_manual(values = mycols_founder)+
  ylim(0,1000)+
  theme_classic()+
  xlab('Subtype')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")

# 1.4 Number of consensus sequences analysed
numseqs_df <- GetNumSeqs(df)

p1.numseq <- ggplot(numseqs_df , aes(sequencing.number))+
  geom_histogram(binwidth=5, (aes(fill = forcats::fct_rev(factor(multiple.founders)), colour = forcats::fct_rev(factor(multiple.founders)))))+
  scale_color_manual(values = mycols_founder)+
  scale_fill_manual(values = mycols_founder)+
  ylim(0,1000)+
  theme_classic()+
  xlab('Number of Consensus Genomes Analysed')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency') +
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity")

# 1.5 Route of transmission (reported.exposure)
names <- c('reported.exposure' , 'multiple.founders' , 'frequency')
exposure_grouped <- cbind.data.frame(df$reported.exposure, df$multiple.founders) %>% split_transission(. ,names)
colnames(exposure_grouped ) <- c('reported.exposure' , 'multiple.founders' )

p1.exposure <- ggplot(exposure_grouped, aes(reported.exposure))+
  geom_bar(aes(fill = forcats::fct_rev(factor(multiple.founders)), colour = forcats::fct_rev(factor(multiple.founders))))+
  scale_color_manual(values = mycols_founder)+
  scale_fill_manual(values = mycols_founder)+
  ylim(0,1000)+
  theme_classic()+
  xlab('Reported Exposure')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency')+
  labs(fill = "Founder Multiplicity", colour = "Founder Multiplicity") 

# Combine basic covariate plots into figure (with labels and axis)
legend <- get_legend(p1.seropos + guides(color = guide_legend(nrow = 2)) +
                       theme(legend.position = "bottom"))


tile1.top <- plot_grid(p1.exposure+ theme(legend.position="none"),
                       p1.method+ theme(legend.position="none"),
                       p1.subtype+ theme(legend.position="none"), 
                       p1.seropos+ theme(legend.position="none"),
                       ncol = 2 , nrow = 2,align = "hv", axis = "bt" , labels = "AUTO") 

tile1.bottom <- plot_grid(p1.numseq+ theme(legend.position="none"),
                          legend, ncol = 2 , nrow = 1, labels = c("E"),rel_widths =c(1, 1, .3)) 


figure1 <- cowplot::plot_grid(tile1.top, 
                              tile1.bottom, 
                              align = "hv", axis = "bt", nrow= 2, ncol = 1, rel_heights =  c(2, 1))

# Print to file


###################################################################################################
###################################################################################################
#bar plot for exposure
#Risk
names <- c('reported.exposure' , 'sub.exposure' , 'frequency')
exposures_df <- stacked_categories(df$reported.exposure, names)


p6 <- ggplot(exposures_df, aes(x = reported.exposure , y = frequency))+
  geom_bar(stat = 'identity' , aes(fill = sub.exposure, colour = sub.exposure), position = 'stack')+
  scale_color_manual(values = mycols_method)+
  scale_fill_manual(values = mycols_method)+
  theme_classic()+
  xlab('Exposure')+
  theme( axis.text.x=element_text(angle=45, hjust=1))+
  ylab('Frequency') + labs(fill = "Reported Exposure", colour = "Reported Exposure")

p6


#Year of Publication

#Minimum numcitqtber of founders plot





###others
df <- read.csv("yop.csv")
dates <- as.Date(ISOdate(df$Year.of.Publication, 1, 1))
df.yrs <- as.data.frame(dates)
start <- as.Date('1991-6-30')
end <- as.Date('2020-6-30')
colnames(df.yrs) <- 'Year.of.Publication'
p1 <- ggplot(df.yrs, aes(Year.of.Publication ))+
  geom_bar(fill = '#000066')+
  geom_vline(xintercept = as.numeric(dates[58]),linetype=2, colour="red" )+
  scale_color_brewer(palette = 'RdBu')+
  scale_fill_brewer(palette = 'RdBu')+
  theme_classic()+
  ylab('')+
  xlab('Year of Publication')+scale_x_date(date_labels = "%Y", date_breaks = '1 year' , limits = c(start,end))+
  theme(axis.text.x=element_text(angle=45 , hjust = 1))

p1


df <- read.csv("subtype.csv")
dates <- as.Date(ISOdate(df$Year.of.Publication, 1, 1))
df.yrs <- as.data.frame(dates)
start <- as.Date('1991-6-30')
end <- as.Date('2020-6-30')
colnames(df.yrs) <- 'Year.of.Publication'
p1 <- ggplot(df.yrs, aes(Year.of.Publication ))+
  geom_bar(fill = '#000066')+
  geom_vline(xintercept = as.numeric(dates[58]),linetype=2, colour="red" )+
  scale_color_npg()+
  scale_fill_npg()+
  theme_classic()+
  ylab('')+
  xlab('Year of Publication')+scale_x_date(date_labels = "%Y", date_breaks = '1 year' , limits = c(start,end))+
  theme(axis.text.x=element_text(angle=45 , hjust = 1))

p1




