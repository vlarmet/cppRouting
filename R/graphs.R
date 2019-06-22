#' Construct graph
#' 
#' @param df  A data.frame or matrix containing 3 columns: from, to, cost. See details.
#' @param directed logical. If FALSE, then all edges are duplicated by inverting 'from' and 'to' nodes.
#' @param coords Optional. A data.frame or matrix containing all nodes coordinates. Columns order should be 'node_ID', 'X', 'Y'.
#' @return List
#' @details 'from' and 'to' are character or numeric vector containing nodes IDs. 
#' 'cost' is a non-negative numeric vector describing the cost (e.g time, distance) between each 'from' and 'to' nodes.
#' coords should not be angles (e.g latitude and longitude), but expressed in a projection system. 
#' 
#' @examples 
#' #Data describing edges of the graph
#' edges<-data.frame(from_vertex=c(0,0,1,1,2,2,3,4,4), 
#'                   to_vertex=c(1,3,2,4,4,5,1,3,5), 
#'                   cost=c(9,2,11,3,5,12,4,1,6))
#'                   
#' #Construct directed and undirected graph 
#' directed_graph<-makegraph(edges,directed=TRUE)
#' non_directed<-makegraph(edges,directed=FALSE)
#' 
#' #Visualizing directed and undirected graphs
#' if(requireNamespace("igraph",quietly = TRUE)){
#'   plot(igraph::graph_from_data_frame(edges))
#'   plot(igraph::graph_from_data_frame(edges,directed=FALSE))
#' } 
#' 
#' #Coordinates of each nodes
#' coord<-data.frame(node=c(0,1,2,3,4,5),X=c(2,2,2,0,0,0),Y=c(0,2,2,0,2,4))
#' 
#' #Construct graph with coordinates
#' directed_graph2<-makegraph(edges, directed=TRUE, coords=coord)
#' 
#' 
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
  if (any(df[,3]<0)) stop("Negative cost is not allowed")
  if (sum(df[,1]==df[,2])>0) df<-df[df[,1]!=df[,2],]
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
