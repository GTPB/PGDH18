
####Histogram####


Hist_ibd<-function(obs.file,files,marker,file.path, project,path_to_save=NULL){
  
  require(reshape2)
  require(ggplot2)
  require(plyr)
  
  
  grab(files,file.path,type="ibd_raw",marker) #It will take the file "genind" from the scenario of interest.
  
  grab(obs.file,file.path,type="ibd_raw",marker)
  
  for(t in 1:length(marker)){
  ibd.obs<-get(paste0(obs.file,"_",marker[t],".ibd_raw"))
  ibd.obs$marker<-marker[t]
  
  if(unique(ibd.obs$marker)=="MT"){
    b="obs. mtDNA"
  }
  if(unique(ibd.obs$marker)=="Y"){
    b="obs. Y chromosome"
  }
  
  ibd.obs$marker<- factor(ibd.obs$marker, labels = b)
  assign(paste0("ibd_obs_line_",obs.file,"_",marker[t]), ibd.obs) #It assigns the dataframe to a unique name which identifies the scenario.
  }
  
  
  for(i in 1:length(files)){
    
    for(t in 1:length(marker)){
      
      ibd<-get(paste0(files[i],"_",marker[t],".ibd_raw"))
      
      
      rate<-grep(c("ibd.rate"), names(ibd))
      inter<-grep(c("ibd.inter"), names(ibd))
      mantel<-grep(c("mantel"), names(ibd))
      
      ibd.dat<-ibd[,-mantel]
      
      ibd.info<-names(ibd.dat)[c(rate,inter)]#It selects the IBD slope and intercept
      
      
      #ibd.dat<-melt(ibd.info)
      ibd.df<-melt(ibd.dat,measure.vars=ibd.info)
      ibd.df$stat<-c(rep("Slope",length(grep("ibd.rate",ibd.df$variable))),rep("Intercept",length(grep("ibd.inter",ibd.df$variable))))
      
      
      
      ibd.df$scen<-files[i] #It adds a column with an abbreviation for each scenario. It will be useful for the plotting part.
      #ibd.obs.df$scen<-files[i]
      
      ibd.df$marker<-marker[t]
     
      cux<-strsplit(files[i], split= "\\_")[[1]]
      
      ibd.df$migr<-cux[grep("msr",cux)]
      
      
      if(unique(ibd.df$migr)=="msr01"){
        a="Female philopatry"
      }
      if(unique(ibd.df$migr)=="msr05"){
        a="Bilocality"
      }
      if(unique(ibd.df$migr)=="msr09"){
        a="Male philopatry"
      }
      
      ibd.df$migr <- factor(ibd.df$migr, labels = a)
      
      if(unique(ibd.df$marker)=="MT"){
        b="mtDNA"
      }
      if(unique(ibd.df$marker)=="Y"){
        b="Y chromosome"
      }
      
     
      ibd.df$marker<- factor(ibd.df$migr, labels = b)
      
      
      assign(paste0("ibd_histo_",files[i],"_",marker[t]), ibd.df) #It assigns the dataframe to a unique name which identifies the scenario.
    }
    
  }
  
  
  
  files.frame<-ls(pattern="ibd_histo_") #It lists the objects with the specified pattern.
  ibd.hist<-do.call(rbind,mget(files.frame)) #It merges the different files with this pattern.
  rm(list=ls(pattern="ibd_histo_"))
  
  files.obs<-ls(pattern="ibd_obs_line_") #It lists the objects with the specified pattern.
  ibd.obs.hist<-do.call(rbind,mget(files.obs)) #It merges the different files with this pattern.
  rm(list=ls(pattern="ibd_obs_line_"))
  
  obs.hist<-ibd.obs.hist[,-c(grep("mantel", names(ibd.obs.hist)))]
  
  obs.df.hist<-melt(obs.hist,measure.vars=c("ibd.rate.1","ibd.inter.1"))
  obs.df.hist$stat<-c(rep("Slope",length(grep("ibd.rate",obs.df.hist$variable))),rep("Intercept",length(grep("ibd.inter",obs.df.hist$variable))))
  
  plot1<-ggplot(ibd.hist, aes(value, fill = marker)) + geom_density(alpha = 0.2,aes(linetype=stat))+
    facet_grid(stat ~ migr,scales="free")+labs(list(fill="Marker",linetype=""))+
    geom_vline(data=ddply(obs.df.hist, stat ~ marker),mapping=aes(xintercept=value,linetype=stat,color=marker))+
    theme_bw()
  
  if(is.character(path_to_save)==TRUE){
    ggsave(plot = plot1,filename = paste0(path_to_save,"/",project,".png"),width = 300, height = 200 ,device = "png", units="mm")
  }else
  {
    print(plot1)
  }
  
}

