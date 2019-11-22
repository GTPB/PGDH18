##marker: codom; mtDNA; Y_chrom.
##Still to be complemented with the function to compute fst with mtDNA###


ibd.intime<-function(files, file.path,marker, save){
  
  require("ade4")
  
  grab(files,file.path,type="fst_raw",marker) #It will take the file "genind" from the scenario of interest.
  
  for(j in 1:length(files)){
    
    fst_raw<-get(paste0(files[j],"_",marker[j],".fst_raw"))
    
    n<-unique(fst_raw$gen)#Get the list of sampled generations
    sim<-length(grep("Sim", names(fst_raw))) #number of simulations present in the output folder
    
    for(p in 1:sim){
      mod<-fst_raw[,which(names(fst_raw)==paste0("Sim",p))]
      r<-as.numeric(gsub("NaN",0,mod))
      fst_raw[[paste0("Sim",p)]]<-r
      
    }
    
    
    ibd.tot.inter<-NULL
    ibd.tot.slope<-NULL
    mantel.tot<-NULL
    
    for(i in 1:length(n)){
      ibd.inter<-NULL
      ibd.slope<-NULL
      mantel.r<-NULL
      mantel.pval<-NULL
      
      for(k in 1:sim){
        data.gen<-fst_raw[which(fst_raw[,"gen"]==n[i]),] #Subset values for a specific timepoint
        linear<-lm(unlist(data.gen[,paste0("Sim",k)]) ~ log(data.gen[,"geo.dist"])) #Compute ibd at a specific timepoint
        ibd.slope<-c(ibd.slope, linear$coefficients[2])
        ibd.inter<-c(ibd.inter, linear$coefficients[1])
        mantel.gen<-data.gen
        nams<-with(mantel.gen, unique(c(as.character(popA), as.character(popB))))
        nams.geo<-with(mantel.gen, unique(c(as.character(popA), as.character(popB))))
        tut<-mantel.gen[,grep(paste0("Sim",k), names(mantel.gen))]
        dist.gen<-structure(tut, Size = length(nams), Labels = nams,Diag = FALSE,Upper = FALSE,method = "user",class = "dist")
        dist.geo<-with(mantel.gen, structure(geo.dist, Size = length(nams.geo), Labels = nams.geo,Diag = FALSE,Upper = FALSE,method = "user",class = "dist"))
        manty<-mantel.randtest(dist.gen, dist.geo, nrepet = 9999)
        mantel.pval<-c(mantel.pval, manty$pvalue) #record Mantel r value and p-value
        mantel.r<-c(mantel.r,manty$obs)
        mantel.sim<-c(mantel.r, mantel.pval)
      }
      
      ibd.tot.slope<-rbind(ibd.tot.slope,ibd.slope)
      ibd.tot.inter<-rbind(ibd.tot.inter,ibd.inter)
      mantel.tot<-rbind(mantel.tot, mantel.sim)
      
      
    }
    
    
    ibd.data.frame<-cbind(n,ibd.tot.slope, ibd.tot.inter,mantel.tot)
    colnames(ibd.data.frame)<-c("gen",paste0("ibd.rate.", seq(1, sim, 1)), paste0("ibd.inter.", seq(1, sim, 1)),paste0("mantel.r.", seq(1, sim, 1)), paste0("mantel.pval.", seq(1, sim, 1)))
    
    if(isTRUE(save)){
      
      write.csv(as.matrix(ibd.data.frame), file=paste0(file.path,"/",files[j],"/",files[j],"/",files[j],"_ibd_raw_",marker[j],".csv"), row.names=FALSE)
      
      
    }else{
      return(ibd.data.frame)
    }
    
    
  }
}


