
####Histogram####





Hist_ibd<-function(files,marker,file.path, project,path_to_save=NULL){

 require(reshape2)
 require(ggplot2)
 require(plyr)
  
  
  grab(files,file.path,type="ibd_raw",marker) #It will take the file "genind" from the scenario of interest.
  


  for(i in 1:length(files)){
    
    ibd<-get(paste0(files[i],"_",marker[i],".ibd_raw"))

    
    rate<-grep(c("ibd.rate"), names(ibd))
    inter<-grep(c("ibd.inter"), names(ibd))
    
    ibd.info<-ibd[,c(rate,inter)]#It selects the IBD slope and intercept
    
    
    ibd.dat<-melt(ibd.info)
    
    ibd.dat$scen<-files[i] #It adds a column with an abbreviation for each scenario. It will be useful for the plotting part.
    ibd.dat$marker<-marker[i]
    cux<-strsplit(files[i], split= "\\_")[[1]]
    ibd.dat$migr<-cux[grep("msr",cux)]
    
    ibd.dat$obs[which(ibd.dat$variable=="ibd.rate.1")]<-0.06117475
    ibd.dat$obs[which(ibd.dat$variable=="ibd.inter.1")]<-0.4233067
    ibd.dat$variable <- factor(ibd.dat$variable, labels = c("Slope", "Intercept"))
    
    if(unique(ibd.dat$migr)=="msr01"){
      a="Female philopatry"
    }
    if(unique(ibd.dat$migr)=="msr05"){
      a="Bilocality"
    }
    if(unique(ibd.dat$migr)=="msr09"){
      a="Male philopatry"
    }
    
    ibd.dat$migr <- factor(ibd.dat$migr, labels = a)
    
    if(unique(ibd.dat$marker)=="MT"){
      b="mtDNA"
    }
    if(unique(ibd.dat$marker)=="Y"){
      b="Y chromosome"
    }
    
    ibd.dat$marker<- factor(ibd.dat$migr, labels = b)
    
    
    assign(paste0("ibd_histo_",files[i],"_",marker[i]), ibd.dat) #It assigns the dataframe to a unique name which identifies the scenario.
    
  }
  
  
  
  
  files.frame<-ls(pattern="ibd_histo_") #It lists the objects with the specified pattern.
  ibd.hist<-do.call(rbind,mget(files.frame)) #It merges the different files with this pattern.
  rm(list=ls(pattern="ibd_histo_"))
  
  

  
  plot1<-ggplot(ibd.hist, aes(value, fill = marker)) + geom_density(alpha = 0.2,aes(linetype=variable))+
    facet_grid(variable ~ migr,scales="free")+labs(list(fill="Marker",linetype=""))+
    geom_vline(data=ddply(ibd.hist, variable ~ migr), 
                 mapping=aes(xintercept=obs,linetype=variable), color="red")
   
  
  
  
  if(is.character(path_to_save)==TRUE){
    pdf(paste0(path_to_save,"/",project,".pdf"))
    par(mfrow=c(1,1))
    print(plot1)
    dev.off()
  }else
  {
    print(plot1)
  }
  
}

