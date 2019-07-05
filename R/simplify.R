cpp_simplify<-function(Graph,keep=NULL,new_edges=FALSE,rm_loop=TRUE,iterate=FALSE,silent=TRUE){
  
  #Nodes to keep
  to_keep<-rep(0,Graph$nbnode)
  if (!is.null(keep)) {
    to_keep[Graph$dict$ref %in% keep]<-1
  }
  
  
  simp<-Simplify2(Graph$data$from,Graph$data$to,Graph$data$dist,Graph$nbnode,loop=rm_loop,keep = to_keep,dict = Graph$dict$ref)
  
  
  
  if (new_edges==TRUE)edges<-list(simp[[3]])
  else edges<-NULL

  
  
  
  #Because removing nodes can create other nodes to remove
  counter<-1
  while(iterate==TRUE){
    if (counter==1 & silent==FALSE) message(paste("  iteration :",counter,"-",Graph$nbnode-simp[[2]],"nodes removed"))
    
    if (simp[[2]]==Graph$nbnode) break
    count<-simp[[2]]
    
    
    
    rd<-Remove_duplicate(simp[[1]][,1],simp[[1]][,2],simp[[1]][,3],Graph$nbnode)
    simp<-Simplify2(rd[,1],rd[,2],rd[,3],Graph$nbnode,loop=rm_loop,keep = to_keep,dict = Graph$dict$ref)
    
    counter<-counter+1
    
    
    
    
    if (count == simp[[2]]) break
    
    if(silent==FALSE) message(paste("  iteration :",counter,"-",count-simp[[2]],"nodes removed"))
    
    
    
    if (new_edges==TRUE) edges[[length(edges)+1]]<-simp[[3]]
    
  }
  
  rd<-Remove_duplicate(simp[[1]][,1],simp[[1]][,2],simp[[1]][,3],Graph$nbnode)

  
  simp<-rd
  if (nrow(simp)==0) stop("All nodes have been removed")
  
  Nodes=unique(c(simp[,1],simp[,2]))

  
  dict<-Graph$dict[Graph$dict$id %in% Nodes,]
  dict$idnew<-0:(nrow(dict)-1)
  simp[,1]<-dict$idnew[match(simp[,1],dict$id)]
  simp[,2]<-dict$idnew[match(simp[,2],dict$id)]
  simp<-as.data.frame(simp)
  simp[,1]<-as.integer(simp[,1])
  simp[,2]<-as.integer(simp[,2])
  colnames(simp)<-c("from","to","dist")
  if (!is.null(Graph$coords)){
    coords<-Graph$coords
    coords<-coords[match(dict$id,Graph$dict$id),]
  }
  else coords=NULL
  
  dict<-dict[,-2]
  colnames(dict)<-c("ref","id")
  


  
  return(list(graph=list(data=simp,
                         coords=coords,
                         nbnode=length(Nodes),
                         dict=dict),
              new_edges=edges))
  
  }
