###Script to compare the IBD pattern for the observed scenario ("real study case") with that of
###putative scenarios explaining this pattern
###Please modify the paths (file.path and path_to_save) according to where do you locate them.
###marker: "MT" for mitochondrial, "Y" for Y-chromosome.


obs.file<-"obs"
files<-c("m01_k100_msr01_malesux01","m01_k100_msr05_malesux01","m01_k100_msr09_malesux01")
marker<-c("MT","Y")
file.path<-"C:/Users/gsgarlata/Desktop/Progetti/2017/Classes_2/Practicals/New_Sims/output"
project<-"test"
path_to_save<-"C:/Users/gsgarlata/Desktop/Progetti/2017/Classes_2/Practicals/New_Sims/output"


path_to_functions<-"C:/Users/gsgarlata/Desktop/Progetti/2017/Classes_2/Practicals/R_functions_classes_2018"
source(paste0(path_to_functions,"/","Get_Data_Seq.R"))
source(paste0(path_to_functions,"/","Fst_seq.R"))
source(paste0(path_to_functions,"/","grab.R"))
source(paste0(path_to_functions,"/","Hist_ibd.R"))
source(paste0(path_to_functions,"/","ibd.intime.R"))

Hist_ibd(obs.file,files,marker,file.path, project,path_to_save)

