




IBD_data_analysis<-function(files,marker,layer_ID,generations_list,file.path,path_to_functions,project,path_to_save){

  if (suppressPackageStartupMessages(!require("strataG")))
  {install.packages("strataG", dependencies = TRUE)}  
  if (suppressPackageStartupMessages(!require("reshape2")))
  {install.packages("reshape2", dependencies = TRUE)}  
  if (suppressPackageStartupMessages(!require("reshape")))
  {install.packages("reshape", dependencies = TRUE)}  
  if (suppressPackageStartupMessages(!require("ggplot2")))
  {install.packages("ggplot2", dependencies = TRUE)}  

source(paste0(path_to_functions,"/","Get_Data_Seq.R"))
source(paste0(path_to_functions,"/","Fst_seq.R"))
source(paste0(path_to_functions,"/","grab.R"))
source(paste0(path_to_functions,"/","Hist_ibd.R"))
source(paste0(path_to_functions,"/","ibd.intime.R"))

path_to_genind<-file.path

for(i in 1:length(files)){
  
  nameOfSimulation<-files[i]
  
  
  Get_Data_Seq(nameOfSimulation,number_of_simulations=1,generations_list,marker[i], layer_ID,path_to_genind)
  
  
  

  
  
}

Fst_seq(files,marker,file.path,the.generation.list=generations_list,save=TRUE)

ibd.intime(files,file.path,marker,save=TRUE)

Hist_ibd(files,marker,file.path,project,path_to_save)


}


