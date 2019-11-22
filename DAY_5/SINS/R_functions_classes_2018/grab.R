###THIS FUNCTION ALLOWS TO GET DATA (RAW AND PROCESSED) FROM THE SCENARIO FOLDER### 
#marker:codom, Y_chrom, mtDNA# In this argument the user will list the type of marker used in each file specified as inputs.
#Please "marker" has to be list with the same order as the files argument#

grab<-function(files,path,type,marker){
  
  for(i in 1:length(files)){ #Loop over the scenarios of interest
    
    for(t in 1:length(marker)){
    
    if(type=="genind"){
      load(paste0(path,"/",files[i],"/",files[i],"_",marker[t],".Rdata")) #If genind format is specified, it is loaded as R workspace.
      assign(paste0(files[i],"_",marker[t],".",type),all_sim, envir=.GlobalEnv) #It assigns the genind object to the name of the scenario and the type.
    }else{  
      
      stat<-read.csv(paste0(path,"/",files[i],"/",files[i],"_",type,"_",marker[t],".csv")) #If not genind format it is written as .csv file.
      assign(paste0(files[i],"_",marker[t],".",type),stat,envir=.GlobalEnv)} #It assigns the genind object to the name of the scenario and the type.
    
  }
    }
  
}

