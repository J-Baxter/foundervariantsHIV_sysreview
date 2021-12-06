###################################################################################################
###################################################################################################
# Loads dependencies for Inferring the multiplicity of founder variants initiating HIV-1 
# infection: a systematic review and IPD meta-analysis

# Packages should be installed using renv::restore(), which will load the required packages and 
# versions from renv.lock
###################################################################################################
###################################################################################################
require(renv)

status <- renv::status()

if (!status$synchronized){
  warning('dependencies not synchonised')
  sync <- readline(prompt="Do you wish to run renv::restore() and synchronise the local project environment? [Y/n]")
  if(sync == 'Y'){
    renv::restore()
  }else if(sync == 'n'){
    cat('Proceeding to load required packages as currently installed ...')
  }else{
    cat('Please enter a valid response. [Y/n]')
  }
}

required_packages <- dependencies()$Package

for (package in required_packages){
  require(package, character.only = T)
}

installed_packages <- (.packages())

if(all(required_packages %in% installed_packages)){
  cat('All packages installed.')
}else{
  missing_packages <- required_packages[!(required_packages %in% installed_packages)]
  error_message <- sapply(missing_packages, function(x) paste0(x,'\n'))
  stop('The following packages are not installed: \n', error_message)
}

###################################################################################################
# END #
###################################################################################################