get_detour<-function(Graph,from,to,extra=NULL,allcores=FALSE){
  if (length(from)!=length(to)) stop("From and to have not the same length")
  if (is.null(extra)) stop("No extra cost")
  extra<-as.numeric(extra)
  if (extra<=0) stop("extra must be positive")
  if (any(is.na(cbind(from,to)))) stop("NAs are not allowed in origin/destination nodes")
  from<-as.character(from)
  
  to<-as.character(to)
  allnodes<-c(from,to)
  if (sum(allnodes %in% Graph$dict$ref)<length(allnodes)) stop("Some nodes are not in the graph")
  
  from_id<-Graph$dict$id[match(from,Graph$dict$ref)]
  to_id<-Graph$dict$id[match(to,Graph$dict$ref)]
  
  if (allcores==FALSE) res<-Detour(from_id,to_id,Graph$data$from,Graph$data$to,Graph$data$dist,Graph$nbnode,t=extra,Graph$dict$ref)
  else {
    numWorkers <- parallel::detectCores()
    cl <- parallel::makeCluster(numWorkers, type = "PSOCK")
    parallel::clusterEvalQ(cl = cl,library("cppRouting"))
    chunks <- parallel::splitIndices(length(from), ncl = numWorkers)
    mylist<-lapply(chunks,function(x) from_id[x])
    mylist2<-lapply(chunks,function(x) to_id[x])
    
    
    res<-parallel::clusterMap(cl,Detour,dep=mylist,arr=mylist2,
                              MoreArgs = list(gfrom=Graph$data$from,gto=Graph$data$to,gw=Graph$data$dist,NbNodes=Graph$nbnode,t=extra,dict=Graph$dict$ref))
    parallel::stopCluster(cl)
    
    res<-do.call(c,res)
  }
  
  names(res)<-paste0(from,"_",to)
  return(res)
}

