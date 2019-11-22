
Get_Data_Seq<-function(nameOfSimulation,number_of_simulations,generations_list,marker, layer_ID,path_to_genind){
  
  all_sim = vector(mode = "list", length = number_of_simulations)
  

  
  for(simID in 1:number_of_simulations){
    
    setwd(paste0(file.path,"/",nameOfSimulation,"/",nameOfSimulation, "/Adegenet_sim_", simID, "/"))
    
    raw_data_per_sim = vector(mode = "list", length = length(generations_list))
    print(paste("Simulation ",simID," of ",number_of_simulations, sep=""))
    
    for(generationID in 1:length(generations_list)){
      

      currentGeneration = generations_list[generationID]
      
      name_of_file = paste(currentGeneration,layer_ID,sep = "_")
      
      
      raw_data_per_sim[[generationID]] = as.matrix(read.csv(paste(name_of_file,"_seqSNP_",marker,".raw",sep = ""),
                                                            header=FALSE,sep="\t")) 
      #It reads the raw data for MT and Y
      
      
      ind_names = raw_data_per_sim[[generationID]][-1,2]
      
      
      pop_names = raw_data_per_sim[[generationID]][-1,1]
      
      
      first<-cbind(ind_names,pop_names)
      
      ##Identify haplotypes##
      
      ind_seq<-raw_data_per_sim[[generationID]][-1,7] 
      
      haps<-unique(ind_seq)
      
      haps_numb<-sprintf("%02d", 1:length(haps))
      
      haps_seq<-matrix(haps,ncol=1)
      
      row.names(haps_seq)<-rep(paste0("Hap.",haps_numb))
      
      haplo<-NULL
      
      for( j in 1:length(ind_seq)){
        hap_index<-which(haps_seq==ind_seq[j])
        
        hap_id<-names(haps_seq[hap_index,1])
        
        haplo<-c(haplo,hap_id)
      }
      
      final<-as.data.frame(cbind(first,haplo))
      
      colnames(final)<-c("id","pop","gene")
      
      
      xy_coords = NULL
      
      for(n in 1:length(ind_names)){
        xy_coords= rbind(xy_coords,
                         cbind(unlist(strsplit(ind_names[n],split = "_"))[4], #x position
                               unlist(strsplit(ind_names[n],split = "_"))[5])) #y position
      }
      xy_coords = apply(xy_coords, 2, as.numeric)
      
      rownames(xy_coords)<-pop_names
      colnames(xy_coords)<-c("x","y")
      
      
      
      dl.g<-df2gtypes(final, ploidy = 1, sequences = NULL,description=as.character(currentGeneration),other=xy_coords)
      
      raw_data_per_sim[[generationID]] = dl.g
      
    }
    
    all_sim[[simID]] = raw_data_per_sim
    
  }  
  
  
  if(is.character(path_to_genind)==TRUE){
    save(all_sim,file=paste0(path_to_genind,"/",nameOfSimulation,"/",nameOfSimulation,"/",nameOfSimulation,"_",marker,".Rdata"))
  }else{}
  
  
}



