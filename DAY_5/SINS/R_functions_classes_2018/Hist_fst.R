
####Histogram####



marker<-c("MT","Y")

Hist_fst<-function(obs.file,files,marker,file.path, project,path_to_save=NULL){
  
  require(reshape2)
  require(ggplot2)
  require(plyr)
  
  
  grab(files,file.path,type="fst_raw",marker) #It will take the file "genind" from the scenario of interest.
  
  grab(obs.file,file.path,type="fst_raw",marker)
  
  
  for(i in 1:length(files)){
    
    for(t in 1:length(marker)){
    
    fst<-get(paste0(files[i],"_",marker[t],".fst_raw"))
    
    fst.obs<-get(paste0(obs.file,"_",marker[t],".fst_raw"))

    sims<-grep("Sim", names(fst))

    fst.df<-melt(fst,measure.vars=names(fst)[sims])
    
    sims.obs<-grep("Sim", names(fst.obs))
    
    fst.obs.df<-melt(fst.obs,measure.vars=names(fst.obs)[sims.obs])
    
    fst.df$scen<-files[i] #It adds a column with an abbreviation for each scenario. It will be useful for the plotting part.
    fst.obs.df$scen<-files[i]
    
    fst.df$marker<-marker[t]
    fst.obs.df$marker<-marker[t]
    
    cux<-strsplit(files[i], split= "\\_")[[1]]
   
     fst.df$migr<-cux[grep("msr",cux)]
    
    
    if(unique(fst.df$migr)=="msr01"){
      a="Female philopatry"
    }
    if(unique(fst.df$migr)=="msr05"){
      a="Bilocality"
    }
    if(unique(fst.df$migr)=="msr09"){
      a="Male philopatry"
    }
    
    fst.df$migr <- factor(fst.df$migr, labels = a)
    
    if(unique(fst.df$marker)=="MT"){
      b="mtDNA"
    }
    if(unique(fst.df$marker)=="Y"){
      b="Y chromosome"
    }
    
    if(unique(fst.obs.df$marker)=="MT"){
      b="mtDNA"
    }
    if(unique(fst.obs.df$marker)=="Y"){
      b="Y chromosome"
    }
    
    fst.df$marker<- factor(fst.df$migr, labels = b)
    
    fst.obs.df$marker<- factor(fst.obs.df$marker, labels = b)
    
    assign(paste0("fst_histo_",files[i],"_",marker[t]), fst.df) #It assigns the dataframe to a unique name which identifies the scenario.
    assign(paste0("fst_obs_histo_",files[i],"_",marker[t]), fst.obs.df) #It assigns the dataframe to a unique name which identifies the scenario.
    
  }
  
}
  
  
  
  files.frame<-ls(pattern="fst_histo_") #It lists the objects with the specified pattern.
  fst.hist<-do.call(rbind,mget(files.frame)) #It merges the different files with this pattern.
  rm(list=ls(pattern="fst_histo_"))
  
  files.obs<-ls(pattern="fst_obs_histo_") #It lists the objects with the specified pattern.
  fst.obs.hist<-do.call(rbind,mget(files.obs)) #It merges the different files with this pattern.
  rm(list=ls(pattern="fst_obs_histo_"))
  
  
  
  
  plot1<-ggplot(fst.hist, aes(value, fill = marker)) + geom_density(alpha = 0.2)+
    facet_grid(marker~ migr,scales="free")+labs(list(fill="Marker",linetype=""))+
    geom_density(data=fst.obs.hist,aes(value, fill = marker),alpha = 0.2,linetype="dashed")
  
    #geom_vline(data=ddply(ibd.hist, variable ~ migr), 
     #          mapping=aes(xintercept=obs,linetype=variable), color="red")
  
  
  
  
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

