Fst_seq<-function(files,marker,file.path,the.generation.list,save){
  
  grab(files,file.path,type="genind",marker) #It will take the file "genind" from the scenario of interest.
  
  for(j in 1:length(files)){
    
    gtypes.nested.array<-get(paste0(files[j],"_",marker[j],".genind"))
    
    
    geo.dist<-dist(unique(gtypes.nested.array[[1]][[1]]@other), method = "manhattan", diag = FALSE, upper = FALSE)
    
    geo.dist.matrix<-as.matrix(geo.dist)
    geo.dist.partial.df<-melt(geo.dist.matrix)[melt(upper.tri(geo.dist.matrix))$value,]
    names(geo.dist.partial.df)<-c("PopA","PopB","geo.dist")
    
    #### 2 - Create the distance dataframe column by replicating this list of distances ####
    
    ## Create the distance column for a single simulation ##
    geographic.distance.column<-vector("list",length=length(gtypes.nested.array[[1]]))
    geographic.distance.column<-as.data.frame(rep(geo.dist.partial.df$geo.dist,times=length(the.generation.list)))
    names(geographic.distance.column)<-"geo.dist"
    
    ## create population ID column by repeating each of the population IDS a ##
    ## mumber of times = to the number of generation steps ###
    
    popID.single.popA<-as.character(rep(geo.dist.partial.df$PopA,times=length(the.generation.list)))
    popID.single.popB<-as.character(rep(geo.dist.partial.df$PopB,times=length(the.generation.list)))
    population.ID.columns<-data.frame(cbind(popID.single.popA,popID.single.popB))
    names(population.ID.columns)<-c("popA","popB")
    
    
    generation.column<-data.frame(rep(the.generation.list,each=length(geo.dist)))
    names(generation.column)<-"gen"
    
    
    
    Apply_Fst_to_seq<-function(nested_list){
      ##create empty nested list###
      func_list<-vector("list",length=length(nested_list))
      
      for (i in 1:length(nested_list)){
        func_list[[i]]<-vector("list",length=length(nested_list[[i]]))
        
        for (h in 1:length(nested_list[[i]])){
          geno<-nested_list[[i]][[h]]
          phist<-pairwiseTest(geno,nrep=0,stats="phist",quietly = TRUE, model="N")
          
          end<-phist$result[,"PHIst"]
          
          end[end<0]<-0
          
          func_list[[i]][[h]]<-end
          
        }
      }
      return(func_list)
    }
    
    fst<-Apply_Fst_to_seq(gtypes.nested.array)
    
    temp.df<-lapply(fst,as.data.frame)
    
    temp.df.melt<-lapply(temp.df,melt)
    
    fst.columns<-data.frame(matrix(nrow=length(temp.df.melt[[1]]$value),ncol = length(temp.df.melt)))
    
    column.names<-vector("character",length = length(temp.df.melt))
    for(i in 1:length(temp.df.melt)){
      fst.columns[[i]]<-as.list(temp.df.melt[[i]]$value)
      column.names[[i]]<-paste("Sim",i,sep="")
    }
    
    fst.columns[is.na(fst.columns)]<-0
    
    names(fst.columns)<-column.names
    
    fst.data.frame<-cbind(generation.column,population.ID.columns,geographic.distance.column,fst.columns)
    
    rownames(fst.data.frame)<-NULL
    
    
    if(isTRUE(save)){
      write.csv(as.matrix(fst.data.frame), file=paste0(file.path,"/",files[j],"/",files[j],"/",files[j],"_fst_raw_",marker[j],".csv"), row.names=FALSE)
    }else{
      return(fst.data.frame)}
  }
  
  
  
}
