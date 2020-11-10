#expand pooled method data and group according to single/multiple founder thresholds as stipulated in source paper
#Keele et al. PNAS 2008
#102 participants included
#methods include model, distance, phylogenetic and haplotype(ish) categories
#requires csv with column names in the format as follows: participant.var, method1.results, method2.results


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
  
  split <- lapply(col_names, getsplits) %>% unique()
  
  selected <- lapply(split , function(x) df[ , grepl(x, col_names) ])

  names(selected) <- split
  return(selected)
}


#create labelled dataframes with individual methods
labeldfs <- function(listofdfs){
  covar = listofdfs[[1]]
  measures = listofdfs[-1]
  varnames = names(measures)
  
  withparticipant <- lapply(measures, function(x) cbind.data.frame(covar, x))
  
  output <- list()
  for (i in 1:length(withparticipant)){
    
    df <- cbind.data.frame(withparticipant[[i]] , varnames[i])
    names(df)[ncol(df)] <- 'method'
    output[[i]] <- df
  }
  
  names(output) <- varnames
  return(output)
}


#group
classify <- function(df, criteria){
  require(dplyr)
  require(stats)
  
  multiple <- do.call(dplyr::filter_, list(df, criteria[[1]])) %>% cbind(founder.multiplicity = 'multiple') 
  
  single <- do.call(dplyr::filter_, list(df, criteria[[2]]))  %>% cbind(founder.multiplicity = 'single') 

  nomeasure <- do.call(dplyr::filter_, list(df, criteria[[3]]))  %>% cbind(founder.multiplicity = 'NA') 
  
  classified <- rbind.data.frame(single, multiple, nomeasure)
  stopifnot(nrow(df)==nrow(classified))
  
  return(classified)
  }


#assign multiple/single founder classification according to a constant threshold value
assignclassificationfixed <- function(listofdfs, threshold){
  stopifnot(length(threshold)-1 == length(listofdfs))
  
  #Distance
  distance_df <-  listofdfs$distance
  DISTANCE <- threshold[[1]]
  criteria_d <- c(~ distance.meanpercent >= DISTANCE, ~ distance.meanpercent < DISTANCE, ~ is.na(DISTANCE))
  
  distance_classified <- classify(d, criteria_d)
  
  #Poisson
  poisson_df <- listofdfs$poisson
  POISSON <- threshold[[2]]
  criteria_p <- c(~ poisson.GOF >= POISSON, ~ poisson.GOF < POISSON, ~ is.na(poisson.GOF))
  
  poisson_classified <- classify(poisson_df, criteria_p)
  stopifnot(nrow(poisson_classified) == nrow(distance_classified))
  
  output <- list(distance_classified, poisson_classified)

  return(output)
}


stratifypooledmethods <- function(data, thresholds){

  labelled_dfs <- groupbycols(keele_combined) %>%
    labeldfs()
  
  #classification for constant values
  classified_dfs <- assignclassification(labelled_dfs[-2], thresholds) #note exclusion
  
  #classification for keele et al TMRCA (relative to particpant.feibig)
  
  
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


#define thresholds as stipulated in keele et al 2008
THRESHOLDS <- c(0.86,0.05)
names(THRESHOLDS) <- c('DISTANCE','POISSON') 



#BEAST TMRCA threshold relative to fiebig stage (ie not a constant)

#run script
stratifypooledmethods(keele_combined, thresholds = THRESHOLDS)

#END#