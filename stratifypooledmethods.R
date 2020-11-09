#expand pooled method data and group according to single/multiple founder thresholds as stipulated in source paper
#Keele et al. PNAS 2008
#102 participants included
#methods include model, distance, phylogenetic and haplotype(ish) categories

groupbycols <- function(df, split){
  col_names = names(df)
  selected = list()
  
  for (i in 1:length(split)){
    
    selected[[i]] = df[ , grepl(split[i], col_names) ]
    
  }
  names(selected) <- split
  return(selected)
}


assignclassification <- function(){
  
}


stratifypooledmethods <- function(df, strata, thresholds){
  grouped_dfs <- groupbycols(df , strata)
  
}


#import data and set groups
keele_combined <- read_csv("keele_combined.csv")
groups <- c('participant', 'distance', 'beast', 'poisson')
thresholds <- data.frame()

#run script
grouped_dfs <- stratifypooledmethods(keele_combined, groups, thresholds)

#END#

#create labelled dataframes with individual methods
distance_df <- cbind.data.frame(grouped_dfs$participant , grouped_dfs$distance)

phylo_df <- cbind.data.frame(grouped_dfs$participant , grouped_dfs$beast)

model_df <- cbind.data.frame(grouped_dfs$participant , grouped_dfs$poisson)

#group dataframes according to thresholds in keele 2008

