# set paths
# << makes it global variable

set_paths<- function(currently_working_env){
  if (currently_working_env == "remote"){
    ## remote data paths (RDSF)
    ## currently not using remotely
    #data_path_tophits <<- "../../data/GWAS_tophits/"  
    #data_path_gwas <<- "../../data/GWAS_tidy/"
    #results_path <<-  "../../results/"
  
  } else if (currently_working_env == "local"){
    ## local data paths
    local_path <<-"/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project4/" #### CHANGE ME
    data_path <<-           paste0(local_path,"01_Data/") 
    data_path_gwas_raw <<-  paste0(local_path,"01_Data/GWAS_raw/") 
    data_path_gwas <<-      paste0(local_path,"01_Data/GWAS_tidy/") 
    data_path_tophits <<-   paste0(local_path,"01_Data/GWAS_tophits/") 
    results_path <<-        paste0(local_path,"02_Results/")
    
  }
  
  # check that all paths exist and create if they don't
  if(!dir.exists(local_path)){ stop("The main project path does not exist")}
  
  dir_path <- gsub(" ", "\\ ", data_path_gwas_raw, fixed=T) # create path vector escaping spaces, otherwise system call cant process it
  if(!dir.exists(dir_path)){ system(paste("mkdir -p", dir_path))}
  dir_path <- gsub(" ", "\\ ", data_path_gwas, fixed=T) # create path vector escaping spaces, otherwise system call cant process it
  if(!dir.exists(dir_path)){ system(paste("mkdir -p", dir_path))}
  dir_path <- gsub(" ", "\\ ", data_path_tophits, fixed=T) # create path vector escaping spaces, otherwise system call cant process it
  if(!dir.exists(dir_path)){ system(paste("mkdir -p", dir_path))}
  dir_path <- gsub(" ", "\\ ", results_path, fixed=T) # create path vector escaping spaces, otherwise system call cant process it
  if(!dir.exists(dir_path)){ system(paste("mkdir -p", dir_path))}
  
}

