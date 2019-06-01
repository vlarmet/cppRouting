#' Construct graph
#' 
#' @param df  A data.frame or matrix containing 3 columns: from, to, cost. See details.
#' @param directed logical. If FALSE, then all edges are duplicated by inverting 'from' and 'to' nodes.
#' @param coords A data.frame or matrix containing all nodes coordinates.
#' @return List
#' @details 'from' and 'to' are character or numeric vector containing nodes IDs. 
#' 'cost' is a numeric vector describing the cost (e.g time, distance) between each 'from' and 'to' nodes.
#' coords should not be angles (e.g latitude and longitude), but expressed in a projection system. 
#' 

makegraph<-function(df,
                    directed=TRUE,
                    coords=NULL){
  df<-as.data.frame(df)
  if (ncol(df)!=3) stop("Data should have 3 columns")
  
  
  df[,1]<-as.character(df[,1])
  df[,2]<-as.character(df[,2])
  df[,3]<-as.numeric(df[,3])
  colnames(df)<-c("from","to","dist")
  
  if (any(is.na(df))) stop("NAs are not allowed in the graph")
  if (sum(df[,1]==df[,2])>0) stop("'Node to the same node' detected")
  #if (sum(duplicated(df))>0) stop("Duplicated vertices not allowed")
  df<-df[!duplicated(df),]
  
  Nodes=unique(c(df[,1],df[,2]))
  
  if (directed==FALSE){
    df2<-df[,c(2,1,3)]
    colnames(df2)<-colnames(df)
    df<-rbind(df,df2)
  }
  
  if (!is.null(coords)){
    if (ncol(coords)!=3) stop("Coords should have 3 columns")
    
    coords[,1]<-as.character(coords[,1])
    coords[,2]<-as.numeric(coords[,2])
    coords[,3]<-as.numeric(coords[,3])
    
    if (any(is.na(coords))) stop("NAs are not allowed in coordinates")
    if (sum(duplicated(coords[,1]))>0) stop("Nodes should be unique in the coordinates data frame")
    #if (length(Nodes)!=nrow(coords)) stop("Number of rows of coords should be equal to the number of egdes")
    if (sum(Nodes %in% coords[,1])<length(Nodes)) stop("Some nodes are missing in coordinates data")
    coords<-coords[coords[,1] %in% Nodes,]
    
  }
  
  dict<-data.frame(ref=Nodes,id=0:(length(Nodes)-1),stringsAsFactors = F)
  
  df[,1]<-dict[match(df[,1],dict$ref),"id"]
  df[,2]<-dict[match(df[,2],dict$ref),"id"]
  coords<-coords[match(Nodes,coords[,1]),]
  
  
  return(list(data=df,
              coords=coords,
              nbnode=length(Nodes),
              dict=dict))
  
}
