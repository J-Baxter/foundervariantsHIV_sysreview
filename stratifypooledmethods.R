#expand pooled method data and group according to single/multiple founder thresholds as stipulated in source paper
#Keele et al. PNAS 2008
#102 participants included
#methods include model, distance, phylogenetic and haplotype(ish) categories


#define functions required in script

getsplits <- function(cols) {
  splitted <- strsplit(cols, '[.]')
  first_element <- unlist(splitted)[1]
  return(first_element)
  }


#split methods into independent dataframes
groupbycols <- function(df){
  
  stopifnot(is(df , 'data.frame'))
  col_names = names(df)
  
  split <- lapply(col_names, getsplits)
  split <- split[!duplicated(split)] #would prefer to use pipe operator here
  
  selected = list()
  for (i in 1:length(split)){
    
    selected[[i]] = df[ , grepl(split[i], col_names) ]
  }
  
  names(selected) <- split
  return(selected)
}


#create labelled dataframes with individual methods
labeldfs <- function(listofdfs, varnames){
  covar = listofdfs[[1]]
  measures = listofdfs[-1]
  
  withparticipant <- lapply(measures, function(x) cbind.data.frame(covar, x))
  
  output <- list()
  for (i in 1:length(withparticipant)){
    
    df <- cbind.data.frame(withparticipant[[i]] , varnames[i+1])
    names(df)[ncol(df)] <- 'method'
    output[[i]] <- df
  }
  
  names(output) <- varnames[-1]
  return(output)
}

 
assignclassification <- function(grouped , threshold){
  #Distance
  
  
  #Poisson
}


stratifypooledmethods <- function(data, varnames, thresholds){
  grouped_dfs <- groupbycols(data , varnames)
  labelled_dfs <- labeldfs(grouped_dfs, varnames)
  classified_dfs <- assignclassification(labelled_dfs, thresholds)
  
  #write output csv(s) to file
  yymmdd <- format(Sys.Date(), '%Y-%b-%d')
  filenames <- lapply(names(classified_dfs) , function(x) paste0('keele', '_', x , '_', yymmdd, '.csv'))
  names(classified_dfs) <- filenames
  
  for (filename in filenames){
    write.csv(filename, file = filename)
  }
}


#import data and set groups
keele_combined <- read_csv("keele_combined.csv")
groups <- c('participant', 'distance', 'beast', 'poisson')

#define thresholds as stipulated in keele et al 2008
THRESHOLDS <- c(0.86,0.05)
names(THRESHOLDS) <- c('DISTANCE','POISSON') 


#BEAST TMRCA threshold relative to fiebig stage (ie not a constant)

#run script
stratifypooledmethods(keele_combined, groups, thresholds = THRESHOLDS)

#END#