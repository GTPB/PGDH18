
##  files: name of the simulated scenario. Input as many times as the type of simulated 
##         genetic markers. (e.g files<- c("test,"test") if mtDNA and Y chromosome
##         have to analysed).
##
##  marker: defines which type of genetic marker has been simulated  (e.g MT and Y).
##
##  layer_ID: defines teh name of the population.
##
##  generations_list: defines the generations for which you have recorded genetic data.   
##
##  file.path: path where output files from SINS are located.
##
##  path_to_functions: path where the R functions required to analyse
##                     simulated genetic data are located.
##
##  project: defines the name of the saved plot, obtained from the analysed data.
##
##  path_to_save: path where the project plot has to be saved.


files<-
marker<-
layer_ID<-
generations_list<-
file.path<-
path_to_functions<-
project<-
path_to_save<-

source(paste0(path_to_functions,"/","Pipeline_classes.R"))


IBD_data_analysis(files,marker,
                  layer_ID,generations_list,
                  file.path,
                  path_to_functions,
                  project,path_to_save)